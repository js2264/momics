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

## Prepare mouse repo for benchmarks

```sh
mkdir tests_data/for_mouse_bench/
wget https://www.encodeproject.org/files/ENCFF763GCB/@@download/ENCFF763GCB.bigWig -O tests_data/for_mouse_bench/MAFK_ENCFF763GCB.bw
wget https://www.encodeproject.org/files/ENCFF880KKC/@@download/ENCFF880KKC.bigWig -O tests_data/for_mouse_bench/CHD2_ENCFF880KKC.bw
wget https://www.encodeproject.org/files/ENCFF857GJE/@@download/ENCFF857GJE.bigWig -O tests_data/for_mouse_bench/H3K4me3_ENCFF857GJE.bw
wget https://www.encodeproject.org/files/ENCFF531JKE/@@download/ENCFF531JKE.bigWig -O tests_data/for_mouse_bench/H4K4me1_ENCFF531JKE.bw
wget https://www.encodeproject.org/files/ENCFF163SBS/@@download/ENCFF163SBS.bigWig -O tests_data/for_mouse_bench/H3K27ac_ENCFF163SBS.bw
wget https://www.encodeproject.org/files/ENCFF705OWT/@@download/ENCFF705OWT.bigWig -O tests_data/for_mouse_bench/H3K9ac_ENCFF705OWT.bw
wget https://www.encodeproject.org/files/ENCFF097KTK/@@download/ENCFF097KTK.bigWig -O tests_data/for_mouse_bench/H3K36me3_ENCFF097KTK.bw
wget https://www.encodeproject.org/files/ENCFF707HHX/@@download/ENCFF707HHX.bigWig -O tests_data/for_mouse_bench/HCFC1_ENCFF707HHX.bw
wget https://www.encodeproject.org/files/ENCFF643WMY/@@download/ENCFF643WMY.bigWig -O tests_data/for_mouse_bench/ZNF384_ENCFF643WMY.bw
wget https://www.encodeproject.org/files/ENCFF550KLM/@@download/ENCFF550KLM.bigWig -O tests_data/for_mouse_bench/ZC3H11A_ENCFF550KLM.bw
./paper/bench/scripts/remove_chrs.py

momics delete -y mm10.momics
momics create mm10.momics
momics ingest chroms -f tests_data/mm10.chrom.sizes -g S288c mm10.momics
momics ingest seq -f tests_data/mm10.fa mm10.momics
momics ingest tracks -f blacklist=tests_data/mm10_blacklist.bed mm10.momics
momics ingest bulk --threads 18 --folder tests_data/for_mouse_bench/fixed/ mm10.momics
momics consolidate mm10.momics
cp -rf mm10.momics tests_data/
momics delete -y mm10.momics
```
