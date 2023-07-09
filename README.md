<div align="center">
<h1 align="center">MAVEQC</h1>
  <p align="center">A R package of quality control on SGE data</p>
</div>

## Table of Contents
<details open>
<summary><b>[Show or Hide]</b></summary>

1. [Dependencies](#dependencies)
2. [Installation](#installation)
3. [File Format](#file-format)
4. [Import Data](#import-data)
5. [Plasmid QC](#plasmid-qc)
    - [QC 1: Sample QC](#pqc1)
    - [QC 2: Experiment QC](#pqc2)
6. [Screen QC](#screen-qc)
    - [QC 1: Sample QC](#sqc1)
    - [QC 2: Experiment QC](#sqc2)
7. [Others](#others)
    - [Test datasets](#test)
    - [Conda](#conda)

</details>

<!-- Dependencies-->
## Dependencies

```R
install.packages("data.table")
install.packages("Ckmeans.1d.dp")
install.packages("reshape2")
install.packages("ggplot2")
install.packages("gplots")
install.packages("plotly")
install.packages("corrplot")
install.packages("ggcorrplot")
install.packages("see")
install.packages("reactable")

install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("DEGreport")
BiocManager::install("apeglm")
```

<p align="right">(<a href="#top">TOP</a>)</p>

<!-- Installation-->
## Installation

Install from github
```R
install.packages("devtools")

library(devtools)
install_github("wtsi-hgi/MAVEQC")
```
Or

Install from source file
```R
install.packages("/path/of/MAVEQC.tar.gz", type = "source")
```

<p align="right">(<a href="#top">TOP</a>)</p>

<!-- File Format-->
## File Format

### sample sheet -- tsv
| sample_name  | library_independent_count | library_dependent_count | valiant_meta | vep_anno | adapt5 | adapt3 | library_name | library_type|
| - | - | - | - | - | - | - | - | - |
| sample1 | s1.allcounts.tsv.gz | s1.libcounts.tsv.gz | meta.csv.gz | meta_consequences.tsv.gz | CTGACTGGCACCTCTTCCCCCAGGA | CCCCGACCCCTCCCCAGCGTGAATG | libA | screen |
| sample2 | s2.allcounts.tsv.gz | s2.libcounts.tsv.gz | meta.csv.gz | meta_consequences.tsv.gz | CTGACTGGCACCTCTTCCCCCAGGA | CCCCGACCCCTCCCCAGCGTGAATG | libA | screen |
| sample3 | s3.allcounts.tsv.gz | s3.libcounts.tsv.gz | meta.csv.gz | meta_consequences.tsv.gz | CTGACTGGCACCTCTTCCCCCAGGA | CCCCGACCCCTCCCCAGCGTGAATG | libA | screen |

* *please use the same headers in the example*
* *adapt5 and adapt3 are required if you don't provide the ref seq and pam seq*
* *vep_anno, library_name and library_type are not necessary, leave them blank if not available*

### library dependent counts -- tsv or tsv.gz
| ID | NAME | SEQUENCE | COUNT | UNIQUE | SAMPLE |
| - | - | - | - | - | - |
| id1 | name1 | ACTTTTCT | 32 | 1 | sample1 | 
| id2 | name2 | ATCTTTCT | 132 | 0 | sample1 | 
| id3 | name3 | ATTCTTCT | 2 | 1 | sample1 | 

* *please use the same headers in the example*
* *please make sure library dependent sequences match with valiant meta file, and sorted by the position in the meta file*

### library independent counts -- tsv or tsv.gz
| SEQUENCE | LENGTH | COUNT |
| - | - | - |
| ACTTTTCT | 276 | 32 | 
| ATCTTTCT | 275 | 132 | 
| ATTCTTCT | 275 | 2 | 

* *please use the same headers in the example*

### valiant meta file
Please use the valiant output file

### vep annotation file
Please use one to one mapping file


<p align="right">(<a href="#top">TOP</a>)</p>

<!-- Import Data -->
## Import Data

### Import from a directory:
All the files are in the same directory including library dependent counts, library independent counts, valiant meta csv, vep annotation and the sample sheet.

#### Load dependencies if required
```
library(data.table)
library(Ckmeans.1d.dp)
library(reshape2)
library(ggplot2)
library(gplots)
library(plotly)
library(corrplot)
library(ggcorrplot)
library(ggbeeswarm)
library(DESeq2)
library(DEGreport)
library(apeglm)
library(see)
library(reactable)
```

#### Import files
```R
library(MAVEQC)

sge_objs <- import_sge_files("/path/to/input/directory", "sample_sheet.tsv")
samqc <- create_sampleqc_object(sge_objs)
```

<p align="right">(<a href="#top">TOP</a>)</p>

<!-- Plasmid QC -->
## Plasmid QC

<a id="pqc1"></a>
### QC 1: Sample QC
```R
samqc <- run_sample_qc(samqc, "plasmid")

qcplot_readlens(samqc, "/path/to/out/plots")
qcplot_stats_total(samqc, "plasmid", "/path/to/out/plots")
qcplot_stats_accepted(samqc, "/path/to/out/plots")
qcplot_position(samqc, "plasmid", "/path/to/out/plots")
```

<p align="right">(<a href="#top">TOP</a>)</p>

<a id="pqc2"></a>
### QC 2: Experimental QC

<p align="right">(<a href="#top">TOP</a>)</p>

<!-- Screen QC -->
## Screen QC

<a id="sqc1"></a>
### QC 1: Sample QC
Reference samples must be assigned. You can use ```select_objects()``` to get this done. ```c(2,5,8)``` is used to point the positions of reference samples in your sample sheet, like 2nd, 5th, 8th line here, or you can use the sample names instead, like ```c("hgsm3_d4_r1","hgsm3_d4_r2","hgsm3_d4_r3")```

```R
samples <- c(2,5,8)
samqc@samples_ref <- select_objects(sge_objs, samples)
samqc <- run_sample_qc(samqc, "screen")

qcplot_readlens(samqc, "/path/to/out/plots")
qcplot_stats_total(samqc, "screen", "/path/to/out/plots")
qcplot_stats_accepted(samqc, "/path/to/out/plots")
qcplot_position(samqc, "screen", "/path/to/out/plots")
qcplot_position_anno(samqc, c("hgsm3_d4_r1", "hgsm3_d4_r2", "hgsm3_d4_r3"), "lof", "/path/to/out/plots")
```

#### coldata example:
| sample_name | replicate | condition |
| - | - | - |
| hgsm3_d4_r1 | R1 | D4 |
| hgsm3_d7_r1 | R1 | D7 |
| hgsm3_d15_r1 | R1 | D15 |
| hgsm3_d4_r2 | R2 | D4 |
| hgsm3_d7_r2 | R2 | D7 |
| hgsm3_d15_r2 | R2 | D15 |
| hgsm3_d4_r3 | R3 | D4 |
| hgsm3_d7_r3 | R3 | D7 |
| hgsm3_d15_r3 | R3 | D15 |


```R
coldata <- read.table("sample_coldata.tsv", header = T, row.names = 1)
expqc <- create_experimentqc_object(samqc, coldata, "D4")
expqc <- run_experiment_qc(expqc)

qcplot_dist_samples(expqc, plotdir = "/path/to/out/plots")
qcplot_pca_samples(expqc, ntop = 500, plotdir = "/path/to/out/plots")
```

<p align="right">(<a href="#top">TOP</a>)</p>

<a id="sqc2"></a>
### QC 2: Experimental QC

```R
qcplot_deseq_fc(expqc, plotdir = "/path/to/out/plots")
```

<p align="right">(<a href="#top">TOP</a>)</p>

<!-- Others -->
## Others

<a id="test"></a>
### Test datasets: 
Test datasets are in the ```test``` directory in the repo. ```plasmid``` and ```screen``` are saperated by analysis. 

<p align="right">(<a href="#top">TOP</a>)</p>

<a id="conda"></a>
### Conda: 
When installing DESeq2, it may have error (Rlog1) on Mac M1 chip. Try cmd below to fix it.

```
export PKG_CPPFLAGS="-DHAVE_WORKING_LOG1P"
```

<p align="right">(<a href="#top">TOP</a>)</p>