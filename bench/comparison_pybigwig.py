from time import sleep
from collections import defaultdict
import pyBigWig
import concurrent.futures
import gc
import multiprocessing
import csv
import math
import subprocess
import pybedtools
import time
import momics
from momics import utils
from momics import multirangequery
from momics.logging import logger
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
            intervals_chrom[
                i * len(intervals_chrom) // ncpus : (i + 1) * len(intervals_chrom) // ncpus
            ]
            for i in range(ncpus)
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

    if os.path.exists("perfs_benchmark_f-nranges.csv"):
        os.remove("perfs_benchmark_f-nranges.csv")

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
            with open("perfs_benchmark_f-nranges.csv", "a", newline="") as f:
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
            with open("perfs_benchmark_f-nranges.csv", "a", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["pybigwig", n_files, n_ranges, w, cpu, end0 - start0, max(mem)])

        if not param.get("skip_pbw") and cpu > 1:
            start0 = time.time()
            mem = memory_usage(lambda: pybigwig_parallel(bws, bed.to_dataframe(), cpu))
            end0 = time.time()
            logger.info(
                f"ELAPSED TIME :: pybigwig_parallel :: {n_files} files, range width {w}, ranges # {n_ranges}, {cpu} threads :: {end0 - start0}"
            )
            with open("perfs_benchmark_f-nranges.csv", "a", newline="") as f:
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


def benchmark_momics_memory_overhead():

    if os.path.exists("perfs_benchmark_f-overhead.csv"):
        os.remove("perfs_benchmark_f-overhead.csv")

    bw_path = "mnase.bw"

    params = [
        {
            "files_n": 10,
            "bin_n": 50000,
            "tile_w": tw,
            "bin_w": w,
            "cpu_n": 16,
        }
        for w in [32, 256, 1024, 2048, 8192, 32768]
        for tw in [50000, 1000, 10000]
    ]
    params = [
        {
            "files_n": 10,
            "bin_n": 500000,
            "tile_w": 50000,
            "bin_w": 32,
            "cpu_n": 16,
        },
        {
            "files_n": 10,
            "bin_n": 500000,
            "tile_w": 1000,
            "bin_w": 32,
            "cpu_n": 16,
        },
        {
            "files_n": 10,
            "bin_n": 500000,
            "tile_w": 10000,
            "bin_w": 32,
            "cpu_n": 16,
        },
        {
            "files_n": 10,
            "bin_n": 1000000,
            "tile_w": 50000,
            "bin_w": 256,
            "cpu_n": 16,
        },
        {
            "files_n": 10,
            "bin_n": 1000000,
            "tile_w": 1000,
            "bin_w": 256,
            "cpu_n": 16,
        },
        {
            "files_n": 10,
            "bin_n": 1000000,
            "tile_w": 10000,
            "bin_w": 256,
            "cpu_n": 16,
        },
        {
            "files_n": 10,
            "bin_n": 250000,
            "tile_w": 50000,
            "bin_w": 1024,
            "cpu_n": 16,
        },
        {
            "files_n": 10,
            "bin_n": 250000,
            "tile_w": 1000,
            "bin_w": 1024,
            "cpu_n": 16,
        },
        {
            "files_n": 10,
            "bin_n": 250000,
            "tile_w": 10000,
            "bin_w": 1024,
            "cpu_n": 16,
        },
        {
            "files_n": 10,
            "bin_n": 32768,
            "tile_w": 50000,
            "bin_w": 8192,
            "cpu_n": 16,
        },
        {
            "files_n": 10,
            "bin_n": 32768,
            "tile_w": 1000,
            "bin_w": 8192,
            "cpu_n": 16,
        },
        {
            "files_n": 10,
            "bin_n": 32768,
            "tile_w": 10000,
            "bin_w": 8192,
            "cpu_n": 16,
        },
        {
            "files_n": 10,
            "bin_n": 8192,
            "tile_w": 50000,
            "bin_w": 32768,
            "cpu_n": 16,
        },
        {
            "files_n": 10,
            "bin_n": 8192,
            "tile_w": 1000,
            "bin_w": 32768,
            "cpu_n": 16,
        },
        {
            "files_n": 10,
            "bin_n": 8192,
            "tile_w": 10000,
            "bin_w": 32768,
            "cpu_n": 16,
        },
    ]
    # 500000 ranges of 2000bp each, for a genome of 12Mb, that represent 1 range every 24bp
    # 10 files, 500,000 ranges of 2,000 bp each = ~ 40Gb RAM
    # 10 files, 100,000 ranges of 8,000 bp each = ~ 32Gb RAM

    for param in params:

        n_files = param["files_n"]
        n_ranges = param["bin_n"]
        tw = param["tile_w"]
        w = param["bin_w"]
        cpu = param["cpu_n"]

        # Create repo
        x = momics.Momics(path=f"test_{n_files}-tracks_overhead.momics").remove()
        x = momics.Momics(path=f"test_{n_files}-tracks_overhead.momics")
        x.add_chroms(utils.get_chr_lengths(bw_path), genome_version="S288c")

        # Add tracks
        bws = {}
        for i in range(n_files):
            bws[f"CH{i}"] = bw_path
        x.add_tracks(bws, threads=18, tile=tw)

        # Prepare bins file for query
        bed = x.bins(w, math.ceil(12000000 / n_ranges), cut_last_bin_out=True)
        print(tw, w, n_ranges, len(bed))

        ## Run momics
        q = multirangequery.MultiRangeQuery(x, bed)
        start = time.time()
        mem = memory_usage(lambda: q.query_tracks(threads=cpu), multiprocess=True)
        end = time.time()
        logger.info(
            f"ELAPSED TIME :: momics :: {n_files} files, tile indexing {tw}, range width {w}, ranges # {n_ranges}, {cpu} threads :: {end - start}"
        )
        with open("perfs_benchmark_f-overhead.csv", "a", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["momics", n_files, n_ranges, tw, w, cpu, end - start, max(mem)])

        del q
        gc.collect()


if __name__ == "__main__":
    # benchmark_momics_pybigwig_nranges()
    benchmark_momics_memory_overhead()
