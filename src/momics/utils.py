import pyBigWig


def get_chr_lengths(bw):
    with pyBigWig.open(bw) as bw:
        a = bw.chroms()
    bw.close()
    return a


def _check_chr_lengths(bw_files, chroms):
    reference_lengths = dict(zip(chroms["chr"], chroms["length"]))
    for file in list(bw_files.values()):
        with pyBigWig.open(file) as bw:
            lengths = bw.chroms()
            if lengths != reference_lengths:
                raise Exception(
                    f"{file} files do not have identical chromomosome lengths."
                )


def _check_track_names(bw_files, tracks):
    labels = set(tracks["label"])
    for element in list(bw_files.keys()):
        if element in labels:
            raise ValueError(
                f"Provided label '{element}' already present in `tracks` table"
            )
