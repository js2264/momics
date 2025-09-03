## Preparing dataset for tests

```sh
momics delete -y test.momics
momics create test.momics
momics ingest chroms -f tests_data/S288c.chrom.sizes -g S288c test.momics
momics ingest seq -f tests_data/S288c.fa test.momics
momics ingest tracks -f ATAC=tests_data/ATAC.bw test.momics
momics ingest tracks -f SCC1=tests_data/SCC1.bw test.momics
momics ingest features -f blacklist=tests_data/S288c_blacklist.bed test.momics
cp -rf test.momics tests_data/
momics delete -y test.momics
```

## Preparing repo for DL

```sh
momics delete -y S288c_MTL.momics
momics create S288c_MTL.momics
momics ingest bulk --folder tests_data/for_multitask_training S288c_MTL.momics
cp -rf S288c_MTL.momics tests_data/
momics delete -y S288c_MTL.momics
```

## Preparing repo for Rossi data

```sh
momics delete -y S288c_Rossi.momics
momics create S288c_Rossi.momics
momics ingest chroms -f tests_data/S288c.chrom.sizes -g S288c S288c_Rossi.momics
momics ingest seq -f tests_data/S288c.fa S288c_Rossi.momics
momics ingest bulk --threads 18 --folder data/bws/fixed/ S288c_Rossi.momics
momics consolidate S288c_Rossi.momics
cp -rf S288c_Rossi.momics tests_data/
momics delete -y S288c_Rossi.momics
```
