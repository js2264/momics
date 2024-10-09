from collections import defaultdict
import numpy as np
import pandas as pd
import pyBigWig
import concurrent.futures
import gc
import tiledb
import multiprocessing
import csv
import math
import pybedtools
import time
import momics
from momics import utils
from momics import multirangequery
from momics.logging import logger
import memory_profiler
from memory_profiler import memory_usage
import os


def _process_batch(bw, intervals):
    results = []
    b = pyBigWig.open(bw)
    for _, row in intervals.iterrows():
        result = b.values(row["chrom"], row["start"], row["end"], numpy=True)
        results.append(result)
    b.close()
    return results


def pybigwig_parallel(bws, intervals, ncpus):

    chroms = list(intervals["chrom"].unique())
    res = []
    for chrom in chroms:
        print(chrom)
        intervals_chrom = intervals[intervals["chrom"] == chrom]
        tasks = [
            intervals_chrom[i * len(intervals_chrom) // ncpus : (i + 1) * len(intervals_chrom) // ncpus] for i in range(ncpus)
        ]
        results_per_track = defaultdict(list)
        for label, bw in bws.items():
            results_per_track[label] = defaultdict(list)
            results_per_track[label][chrom] = []
            with concurrent.futures.ProcessPoolExecutor(max_workers=ncpus) as executor:
                futures = []
                for i in range(ncpus):
                    futures.append(executor.submit(_process_batch, bw, tasks[i]))
                concurrent.futures.wait(futures)

            for future in concurrent.futures.as_completed(futures):
                results_per_track[label][chrom].append(future.result())
        res.append(results_per_track)


def pybigwig(bws, intervals):

    results_per_track = defaultdict(list)
    for label, bw in bws.items():
        results_per_track[label] = _process_batch(bw, intervals)


def benchmark_momics_pybigwig_nranges():

    if os.path.exists("bench/perfs_benchmark_f-nranges.csv"):
        os.remove("bench/perfs_benchmark_f-nranges.csv")

    bw_path = "mnase.bw"

    params = [
        {
            "files_n": 10,
            "bin_n": n,
            "bin_w": 256,
            "cpu_n": cpu,
            "skip_momics": False,
            "skip_bw": False,
            "skip_pbw": False,
        }
        for n in [100, 1000, 10000, 100000, 500000, 1000000, 2500000, 5000000]
        for cpu in [1, 16]
    ]
    # 500000 ranges of 2000bp each, for a genome of 12Mb, that represent 1 range every 24bp
    # 10 files, 500,000 ranges of 2,000 bp each = ~ 40Gb RAM

    for param in params:

        n_files = param["files_n"]
        n_ranges = param["bin_n"]
        w = param["bin_w"]
        cpu = param["cpu_n"]

        # Create repo
        # x = momics.Momics(path=f"test_{n_files}-tracks.momics").remove()
        x = momics.Momics(path=f"test_{n_files}-tracks.momics")
        # x.add_chroms(utils.get_chr_lengths(bw_path), genome_version="S288c")

        # Add tracks
        bws = {}
        for i in range(n_files):
            bws[f"CH{i}"] = bw_path
        # x.add_tracks(bws, threads=18, tile=10000)

        # Prepare bins file for query
        bed = x.bins(w, math.ceil(12000000 / n_ranges), cut_last_bin_out=True)

        ## Run momics
        if not param.get("skip_momics"):
            q = multirangequery.MultiRangeQuery(x, bed)
            start = time.time()
            mem = memory_usage(lambda: q.query_tracks(threads=cpu))
            end = time.time()
            logger.info(
                f"ELAPSED TIME :: momics :: {n_files} files, range width {w}, ranges # {n_ranges}, {cpu} threads :: {end - start}"
            )
            with open("bench/perfs_benchmark_f-nranges.csv", "a", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["momics", n_files, n_ranges, w, cpu, end - start, max(mem)])

            del q
            gc.collect()

        ## Run pybigwig
        if not param.get("skip_bw") and cpu == 1:
            start0 = time.time()
            mem = memory_usage(lambda: pybigwig(bws, bed.to_dataframe()))
            end0 = time.time()
            logger.info(
                f"ELAPSED TIME :: pybigwig :: {n_files} files, range width {w}, ranges # {n_ranges}, {cpu} threads :: {end0 - start0}"
            )
            with open("bench/perfs_benchmark_f-nranges.csv", "a", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["pybigwig", n_files, n_ranges, w, cpu, end0 - start0, max(mem)])

        if not param.get("skip_pbw") and cpu > 1:
            start0 = time.time()
            mem = memory_usage(lambda: pybigwig_parallel(bws, bed.to_dataframe(), cpu))
            end0 = time.time()
            logger.info(
                f"ELAPSED TIME :: pybigwig_parallel :: {n_files} files, range width {w}, ranges # {n_ranges}, {cpu} threads :: {end0 - start0}"
            )
            with open("bench/perfs_benchmark_f-nranges.csv", "a", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(
                    [
                        "pybigwig_parallel",
                        n_files,
                        n_ranges,
                        w,
                        cpu,
                        end0 - start0,
                        max(mem),
                    ]
                )


def _generate_intervals(chroms, N, W):
    chromsizes = {i: j for (i, j) in zip(chroms["chrom"], chroms["length"])}
    intervals = []
    for chrom, size in chromsizes.items():
        step = size // N
        for i in range(N):
            start = i * step + 1
            end = min(start + W, size)
            intervals.append((chrom, start, end))
    intervals = pd.DataFrame(intervals, columns=["chrom", "start", "end"])
    return pybedtools.BedTool.from_dataframe(intervals)


def benchmark_query(tile, chunksize, chrom_length, query_ranges):

    # Schema creation based on tile and chunksize
    dom = tiledb.Domain(
        tiledb.Dim(
            name="position",
            domain=(0, chrom_length),
            dtype=np.int32,
            tile=tile,
            filters=tiledb.FilterList(
                [tiledb.LZ4Filter(), tiledb.ZstdFilter(level=3)],
                chunksize=chunksize,
            ),
        )
    )
    attr = tiledb.Attr(name="placeholder", dtype=np.float32)
    schema = tiledb.ArraySchema(domain=dom, attrs=[attr], sparse=False)

    # Simulate the array
    tdb = "my_benchmark_array"
    tiledb.DenseArray.create(tdb, schema)

    # Write some data
    with tiledb.DenseArray(tdb, mode="w") as A:
        A[:] = np.random.random(chrom_length + 1).astype(np.float32)

    # Time the query
    start_time = time.time()

    # Query either small or large slices
    slices = query_ranges

    # Perform queries
    with tiledb.DenseArray(tdb, mode="r") as A:

        def _run(A, slices):
            A.multi_index[slices]

        mem_usage = memory_profiler.memory_usage(_run(A, slices), interval=0.1)

    # Time and memory reporting
    elapsed_time = time.time() - start_time
    print(f"Tile={tile}, Chunksize={chunksize}: Time: {elapsed_time} s, Max Memory: {max(mem_usage)} MB")


def benchmark_momics_memory_overhead():

    if os.path.exists("bench/perfs_benchmark_f-overhead.csv"):
        os.remove("bench/perfs_benchmark_f-overhead.csv")

    ts0 = 0
    params = [
        # Increasing bin number, decreasing bin width
        {"n_bins": 1000, "bin_width": 40000, "tilesize": 64000},
        {"n_bins": 10000, "bin_width": 4000, "tilesize": 64000},
        {"n_bins": 100000, "bin_width": 400, "tilesize": 64000},
        {"n_bins": 1000000, "bin_width": 40, "tilesize": 64000},
        {"n_bins": 10000000, "bin_width": 4, "tilesize": 64000},
        # Increasing bin number, small bins
        {"n_bins": 100, "bin_width": 40, "tilesize": 64000},
        {"n_bins": 1000, "bin_width": 40, "tilesize": 64000},
        {"n_bins": 10000, "bin_width": 40, "tilesize": 64000},
        {"n_bins": 100000, "bin_width": 40, "tilesize": 64000},
        {"n_bins": 1000000, "bin_width": 40, "tilesize": 64000},
        # Increasing bin number, large bins
        {"n_bins": 100, "bin_width": 4000, "tilesize": 64000},
        {"n_bins": 1000, "bin_width": 4000, "tilesize": 64000},
        {"n_bins": 10000, "bin_width": 4000, "tilesize": 64000},
        {"n_bins": 100000, "bin_width": 4000, "tilesize": 64000},
        {"n_bins": 1000000, "bin_width": 4000, "tilesize": 64000},
        # Same bin number, increasing bins
        {"n_bins": 100000, "bin_width": 4, "tilesize": 64000},
        {"n_bins": 100000, "bin_width": 40, "tilesize": 64000},
        {"n_bins": 100000, "bin_width": 400, "tilesize": 64000},
        {"n_bins": 100000, "bin_width": 4000, "tilesize": 64000},
        {"n_bins": 100000, "bin_width": 40000, "tilesize": 64000},
        # Increasing bin number, decreasing bin width
        {"n_bins": 1000, "bin_width": 40000, "tilesize": 6400000},
        {"n_bins": 10000, "bin_width": 4000, "tilesize": 6400000},
        {"n_bins": 100000, "bin_width": 400, "tilesize": 6400000},
        {"n_bins": 1000000, "bin_width": 40, "tilesize": 6400000},
        {"n_bins": 10000000, "bin_width": 4, "tilesize": 6400000},
        # Increasing bin number, small bins
        {"n_bins": 100, "bin_width": 40, "tilesize": 6400000},
        {"n_bins": 1000, "bin_width": 40, "tilesize": 6400000},
        {"n_bins": 10000, "bin_width": 40, "tilesize": 6400000},
        {"n_bins": 100000, "bin_width": 40, "tilesize": 6400000},
        {"n_bins": 1000000, "bin_width": 40, "tilesize": 6400000},
        # Increasing bin number, large bins
        {"n_bins": 100, "bin_width": 4000, "tilesize": 6400000},
        {"n_bins": 1000, "bin_width": 4000, "tilesize": 6400000},
        {"n_bins": 10000, "bin_width": 4000, "tilesize": 6400000},
        {"n_bins": 100000, "bin_width": 4000, "tilesize": 6400000},
        {"n_bins": 1000000, "bin_width": 4000, "tilesize": 6400000},
        # Same bin number, increasing bins
        {"n_bins": 100000, "bin_width": 4, "tilesize": 6400000},
        {"n_bins": 100000, "bin_width": 40, "tilesize": 6400000},
        {"n_bins": 100000, "bin_width": 400, "tilesize": 6400000},
        {"n_bins": 100000, "bin_width": 4000, "tilesize": 6400000},
        {"n_bins": 100000, "bin_width": 40000, "tilesize": 6400000},
    ]

    for param in params:

        n_bins = param["n_bins"]
        bin_width = param["bin_width"]
        ts = param["tilesize"]
        momics.momics.TILEDB_TILESIZE = ts

        bw_path = "/data/momics/data/mnase_mm10.bw"
        bws = {}
        for i in range(2):
            bws[f"CH{i}"] = bw_path

        if ts != ts0:
            x = momics.Momics(path=f"test_2-tracks_overhead.momics").remove()
            x = momics.Momics(path=f"test_2-tracks_overhead.momics")
            x.add_chroms(utils.get_chr_lengths(bw_path), genome_version="S288c")
            x.add_tracks(bws, threads=23)
            ts0 = ts

        x = momics.Momics(path=f"test_2-tracks_overhead.momics")
        bed = _generate_intervals(x.chroms(), math.ceil(n_bins / len(x.chroms())), bin_width)
        q = multirangequery.MultiRangeQuery(x, bed)

        start = time.time()
        mem = memory_usage(lambda: q.query_tracks(threads=16), interval=0.05, multiprocess=True)
        end = time.time()

        logger.info(
            f"ELAPSED TIME :: momics :: 2 files, tile size {ts}, range width {bin_width}, ranges # {n_bins}, 16 threads :: {end - start}"
        )
        with open("bench/perfs_benchmark_f-overhead.csv", "a", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["momics", 2, ts, n_bins, bin_width, 16, end - start, max(mem)])

        del q
        gc.collect()


def benchmark_momics_coverage_speed_mem_overhead():

    if os.path.exists("bench/perfs_benchmark_f-number_of_file.csv"):
        os.remove("bench/perfs_benchmark_f-number_of_file.csv")

    nf0 = 0
    params = [
        {"gen_cov": 0.000025, "bin_width": 4000, "nfile": 1},
        {"gen_cov": 0.00025, "bin_width": 4000, "nfile": 1},
        {"gen_cov": 0.0025, "bin_width": 4000, "nfile": 1},
        {"gen_cov": 0.025, "bin_width": 4000, "nfile": 1},
        {"gen_cov": 0.25, "bin_width": 4000, "nfile": 1},
        {"gen_cov": 0.000025, "bin_width": 4000, "nfile": 2},
        {"gen_cov": 0.00025, "bin_width": 4000, "nfile": 2},
        {"gen_cov": 0.0025, "bin_width": 4000, "nfile": 2},
        {"gen_cov": 0.025, "bin_width": 4000, "nfile": 2},
        {"gen_cov": 0.25, "bin_width": 4000, "nfile": 2},
        {"gen_cov": 0.000025, "bin_width": 4000, "nfile": 4},
        {"gen_cov": 0.00025, "bin_width": 4000, "nfile": 4},
        {"gen_cov": 0.0025, "bin_width": 4000, "nfile": 4},
        {"gen_cov": 0.025, "bin_width": 4000, "nfile": 4},
        {"gen_cov": 0.25, "bin_width": 4000, "nfile": 4},
        {"gen_cov": 0.000025, "bin_width": 4000, "nfile": 8},
        {"gen_cov": 0.00025, "bin_width": 4000, "nfile": 8},
        {"gen_cov": 0.0025, "bin_width": 4000, "nfile": 8},
        {"gen_cov": 0.025, "bin_width": 4000, "nfile": 8},
        {"gen_cov": 0.25, "bin_width": 4000, "nfile": 8},
        # {"gen_cov": 0.000008, "bin_width": 4000, "nfile": 12},
        # {"gen_cov": 0.00008, "bin_width": 4000, "nfile": 12},
        # {"gen_cov": 0.0008, "bin_width": 4000, "nfile": 12},
        # {"gen_cov": 0.008, "bin_width": 4000, "nfile": 12},
        # {"gen_cov": 0.08, "bin_width": 4000, "nfile": 12},
        # {"gen_cov": 0.000008, "bin_width": 4000, "nfile": 24},
        # {"gen_cov": 0.00008, "bin_width": 4000, "nfile": 24},
        # {"gen_cov": 0.0008, "bin_width": 4000, "nfile": 24},
        # {"gen_cov": 0.008, "bin_width": 4000, "nfile": 24},
        # {"gen_cov": 0.08, "bin_width": 4000, "nfile": 24},
    ]

    for param in params:

        gen_cov = param["gen_cov"]
        bin_width = param["bin_width"]
        nf = param["nfile"]
        GEN_SIZE = 2725537669
        n_bins = math.ceil(gen_cov * GEN_SIZE / bin_width)

        bw_path = "/data/momics/data/mnase_mm10.bw"
        bws = {}
        for i in range(nf):
            bws[f"CH{i}"] = bw_path

        if nf != nf0:
            x = momics.Momics(path=f"test_2-tracks_overhead.momics").remove()
            x = momics.Momics(path=f"test_2-tracks_overhead.momics")
            x.add_chroms(utils.get_chr_lengths(bw_path), genome_version="S288c")
            x.add_tracks(bws, threads=23)
            nf0 = nf

        x = momics.Momics(path=f"test_2-tracks_overhead.momics")
        bed = _generate_intervals(x.chroms(), math.ceil(n_bins / len(x.chroms())), bin_width)
        q = multirangequery.MultiRangeQuery(x, bed)

        start = time.time()
        mem = memory_usage(lambda: q.query_tracks(threads=16), interval=0.05, multiprocess=True)
        end = time.time()

        logger.info(
            f"ELAPSED TIME :: momics :: {nf} files, range width {bin_width}, ranges # {n_bins}, 16 threads :: {end - start}"
        )
        with open("bench/perfs_benchmark_f-number_of_file.csv", "a", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["momics", nf, n_bins, bin_width, gen_cov, end - start, max(mem)])

        del q
        gc.collect()


if __name__ == "__main__":
    multiprocessing.set_start_method("spawn")
    # benchmark_momics_pybigwig_nranges()
    # benchmark_momics_memory_overhead()
    benchmark_momics_coverage_speed_mem_overhead()
