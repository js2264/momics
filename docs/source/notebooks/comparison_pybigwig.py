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

multiprocessing.set_start_method("spawn", True)


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
    for chrom in chroms:
        print(chrom)
        intervals_chrom = intervals[intervals["chrom"] == chrom]
        tasks = [
            intervals_chrom[
                i
                * len(intervals_chrom)
                // ncpus : (i + 1)
                * len(intervals_chrom)
                // ncpus
            ]
            for i in range(ncpus)
        ]
        results_per_track = defaultdict(list)
        for label, bw in bws.items():
            print(label)
            results_per_track[label] = defaultdict(list)
            with concurrent.futures.ProcessPoolExecutor(max_workers=ncpus) as executor:
                futures = []
                for i in range(ncpus):
                    futures.append(executor.submit(_process_batch, bw, tasks[i]))
                concurrent.futures.wait(futures)

            for future in concurrent.futures.as_completed(futures):
                results_per_track[label][chrom] = future.result()


def pybigwig(bws, intervals):

    results_per_track = defaultdict(list)
    for label, bw in bws.items():
        print(label)
        results_per_track[label] = _process_batch(bw, intervals)


def benchmark_momics_pybedtools():

    bw_path = "/home/rsg/Projects/20220309_Christophe_GC-paper/data/WT/ChIP/tracks/CH244/CH244^mapped_S288c^CIJXLY.CPM.bw"

    # params = [
    #     {
    #         "files_n": 10,
    #         "tile_w": 10000,
    #         "bin_n": 100000,
    #         "bin_w": 2000,
    #         "cpu_n": 1,
    #         "skip_pbw": True,
    #     },
    #     {
    #         "files_n": 10,
    #         "tile_w": 10000,
    #         "bin_n": 300000,
    #         "bin_w": 2000,
    #         "cpu_n": 2,
    #         "skip_pbw": True,
    #     },
    #     {"files_n": 10, "tile_w": 10000, "bin_n": 100000, "bin_w": 2000, "cpu_n": 4, "skip_pbw": True},
    #     {"files_n": 10, "tile_w": 10000, "bin_n": 100000, "bin_w": 2000, "cpu_n": 8, "skip_pbw": True},
    #     {"files_n": 10, "tile_w": 10000, "bin_n": 100000, "bin_w": 2000, "cpu_n": 16, "skip_pbw": True},
    #     {"files_n": 10, "tile_w": 10000, "bin_n": 20000, "bin_w": 10, "cpu_n": 1, "skip_pbw": True},
    #     {"files_n": 10, "tile_w": 10000, "bin_n": 20000, "bin_w": 10, "cpu_n": 12, "skip_pbw": True},
    # ]

    params = [
        {
            "files_n": 10,
            "bin_n": n,
            "bin_w": 2000,
            "cpu_n": cpu,
            "skip_momics": False,
            "skip_bw": False,
            "skip_pbw": False,
        }
        for cpu in [16, 4, 1]
        for n in [10000, 100000, 500000]
    ]  # 500000 ranges of 2000bp each, for a genome of 12Mb, that represent 1 range every 24bp

    # params = [
    #     {
    #         "files_n": 10,
    #         "bin_n": nbins,
    #         "bin_w": 2000,
    #         "cpu_n": 1,
    #         "skip_momics": False,
    #         "skip_pbw": True,
    #     }
    #     for nbins in [100000, 250000, 1000000, 2000000]
    # ]

    for param in params:

        n_files = param["files_n"]
        n_ranges = param["bin_n"]
        w = param["bin_w"]
        cpu = param["cpu_n"]

        # Create repo
        x = momics.Momics(path=f"test_{n_files}-tracks.momics").remove()
        x = momics.Momics(path=f"test_{n_files}-tracks.momics")
        x.add_chroms(utils.get_chr_lengths(bw_path), genome_version="S288c")

        # Add tracks
        bws = {}
        for i in range(n_files):
            bws[f"CH{i}"] = bw_path
        x.add_tracks(bws, threads=18, tile=10000)

        # Prepare bins file for query
        bed = x.bins(w, math.ceil(12000000 / n_ranges), cut_last_bin_out=True)

        ## Run momics
        if not param.get("skip_momics"):
            q = multirangequery.MultiRangeQuery(x, bed)
            start = time.time()
            q.query_tracks(threads=cpu)
            end = time.time()
            logger.info(
                f"ELAPSED TIME :: momics :: {n_files} files, range width {w}, ranges # {n_ranges}, {cpu} threads :: {end - start}"
            )
            with open("perfs_benchmark_f-nranges-and-cpu.csv", "a", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["momics", n_files, n_ranges, w, cpu, end - start])

            del q
            gc.collect()

        ## Run pybigwig
        if not param.get("skip_bw") and cpu == 1:
            start0 = time.time()
            q = pybigwig(bws, bed.to_dataframe())
            end0 = time.time()
            logger.info(
                f"ELAPSED TIME :: pybigwig :: {n_files} files, range width {w}, ranges # {n_ranges}, {cpu} threads :: {end0 - start0}"
            )
            with open("perfs_benchmark_f-nranges-and-cpu.csv", "a", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["pybigwig", n_files, n_ranges, w, cpu, end0 - start0])

            del q
            gc.collect()

        if not param.get("skip_pbw") and cpu > 1:
            start0 = time.time()
            q = pybigwig_parallel(bws, bed.to_dataframe(), cpu)
            end0 = time.time()
            logger.info(
                f"ELAPSED TIME :: pybigwig_parallel :: {n_files} files, range width {w}, ranges # {n_ranges}, {cpu} threads :: {end0 - start0}"
            )
            with open("perfs_benchmark_f-nranges-and-cpu.csv", "a", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(
                    ["pybigwig_parallel", n_files, n_ranges, w, cpu, end0 - start0]
                )

            del q
            gc.collect()


if __name__ == "__main__":
    benchmark_momics_pybedtools()
