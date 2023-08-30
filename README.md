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
    - [Report](#pqc-report)
6. [Screen QC](#screen-qc)
    - [QC 1: Sample QC](#sqc1)
    - [QC 2: Experiment QC](#sqc2)
    - [Report](#sqc-report)
7. [Others](#others)
    - [Pandoc](#pandoc)
    - [Test datasets](#test)
    - [Conda](#conda)

</details>

<!-- Dependencies-->
## Dependencies

```R
install.packages("configr")
install.packages("vroom")
install.packages("data.table")
install.packages("Ckmeans.1d.dp")
install.packages("gplots")
install.packages("ggplot2")
install.packages("plotly")
install.packages("ggcorrplot")
install.packages("corrplot")
install.packages("see")
install.packages("ggbeeswarm")
install.packages("reactable")
install.packages("reshape2")
install.packages("htmltools")
install.packages("gtools")

install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("DEGreport")
BiocManager::install("apeglm")
```

Load dependencies if required
```
library(configr)
library(vroom)
library(data.table)
library(Ckmeans.1d.dp)
library(gplots)
library(ggplot2)
library(plotly)
library(ggcorrplot)
library(corrplot)
library(see)
library(ggbeeswarm)
library(reactable)
library(htmltools)
library(reshape2)
library(gtools)

library(DESeq2)
library(DEGreport)
library(apeglm)
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

Install from the compiled source file
```R
install.packages("/path/of/MAVEQC.tar.gz", type = "source")
```

<p align="right">(<a href="#top">TOP</a>)</p>

<!-- File Format-->
## File Format

### sample sheet -- tsv
| sample_name  | library_independent_count | library_dependent_count | valiant_meta | vep_anno | adapt5 | adapt3 | per_r1_adaptor | per_r2_adaptor | library_name | library_type|
| - | - | - | - | - | - | - | - | - | - | - |
| sample1 | s1.allcounts.tsv.gz | s1.libcounts.tsv.gz | meta.csv.gz | meta_consequences.tsv.gz | CTGACTGGCACCTCTTCCCCCAGGA | CCCCGACCCCTCCCCAGCGTGAATG | 0.21 | 0.10 | libA | screen |
| sample2 | s2.allcounts.tsv.gz | s2.libcounts.tsv.gz | meta.csv.gz | meta_consequences.tsv.gz | CTGACTGGCACCTCTTCCCCCAGGA | CCCCGACCCCTCCCCAGCGTGAATG | 0.11 | 0.02 | libA | screen |
| sample3 | s3.allcounts.tsv.gz | s3.libcounts.tsv.gz | meta.csv.gz | meta_consequences.tsv.gz | CTGACTGGCACCTCTTCCCCCAGGA | CCCCGACCCCTCCCCAGCGTGAATG | 0.01 | 0.18 | libA | screen |

* *please use the same headers in the example*
* *adapt5 and adapt3 are required if you don't provide the ref seq and pam seq*
* *vep_anno, library_name and library_type are not necessary, leave them blank if not available*

### library dependent counts -- tsv or tsv.gz
| ID | NAME | SEQUENCE | LENGTH | COUNT | UNIQUE | SAMPLE |
| - | - | - | - | - | - | - |
| id1 | name1 | ACTTTTCT | 276 | 32 | 1 | sample1 | 
| id2 | name2 | ATCTTTCT | 275 | 132 | 0 | sample1 | 
| id3 | name3 | ATTCTTCT | 275 | 2 | 1 | sample1 | 

* *please use the same headers in the example*
* *please make sure library dependent sequences match with valiant meta file*
* *please refer to [pyQUEST](https://github.com/cancerit/pyQUEST#library-dependent-counts)*

### library independent counts -- tsv or tsv.gz
| SEQUENCE | LENGTH | COUNT |
| - | - | - |
| ACTTTTCT | 276 | 32 | 
| ATCTTTCT | 275 | 132 | 
| ATTCTTCT | 275 | 2 | 

* *please use the same headers in the example*
* *please refer to [pyQUEST](https://github.com/cancerit/pyQUEST#library-dependent-counts)*

### valiant meta file
*Please use the VaLiAnT output file, refer to [VaLiAnT](https://github.com/cancerit/VaLiAnt)*

### vep annotation file
*Please use one to one mapping file*


<p align="right">(<a href="#top">TOP</a>)</p>

<!-- Import Data -->
## Import Data

### Import a group of samples from a directory
All the files are in the same directory including library dependent counts, library independent counts, valiant meta csv, vep annotation and the sample sheet.


```R
library(MAVEQC)

sge_objs <- import_sge_files("/path/to/input/directory", "sample_sheet.tsv")
```

<p align="right">(<a href="#top">TOP</a>)</p>

<!-- Plasmid QC -->
## Plasmid QC
> Test datasets are not available now, will add them soon

<a id="pqc1"></a>
### QC 1: Sample QC
```R
output_dir <- "/path/to/output/directory"

samqc <- create_sampleqc_object(sge_objs)
samqc <- run_sample_qc(samqc, "plasmid")

qcplot_samqc_all(samqc, qc_type = "plasmid", plot_dir = output_dir)
qcout_samqc_all(samqc, qc_type = "plasmid", out_dir = output_dir)
```

<p align="right">(<a href="#top">TOP</a>)</p>

<a id="pqc-report"></a>
### Report
This creates a html report concatenating all the results including figures and tables. Please make sure you have generated all the figures and tables, otherwise the report may be incomplete.
```R
create_qc_reports("/path/to/sample/sheet", "plasmid", output_dir)
```

<p align="right">(<a href="#top">TOP</a>)</p>

<!-- Screen QC -->
## Screen QC

<a id="sqc1"></a>
### QC 1: Sample QC
Reference samples must be assigned. You can use ```select_objects()``` to get this done. ```c(2,5,8)``` is used to point the positions of reference samples in your sample sheet, like 2nd, 5th, 8th line here, or you can use the sample names instead, like ```c("hgsm3_d4_r1","hgsm3_d4_r2","hgsm3_d4_r3")```

```R
output_dir <- "/path/to/output/directory"

samqc <- create_sampleqc_object(sge_objs)
ref_samples <- c("hgsm3_d4_r1","hgsm3_d4_r2","hgsm3_d4_r3")
samqc@samples_ref <- select_objects(sge_objs, ref_samples)
samqc <- run_sample_qc(samqc, "screen")

qcplot_samqc_all(samqc, qc_type = "screen", samples = ref_samples, plot_dir = output_dir)
qcout_samqc_all(samqc, qc_type = "screen", out_dir = output_dir)
```

<p align="right">(<a href="#top">TOP</a>)</p>

<a id="sqc2"></a>
### QC 2: Experimental QC
Prepare the coldata for DESeq2 like below.

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
expqc <- create_experimentqc_object(samqc, coldata = coldata, refcond = "D4") 
expqc <- run_experiment_qc(expqc) 

qcplot_expqc_all(expqc, plot_dir = output_dir)
qcout_expqc_all(expqc, out_dir = output_dir)
```

<p align="right">(<a href="#top">TOP</a>)</p>

<a id="sqc-report"></a>
### Report
This creates a html report concatenating all the results including figures and tables. Please make sure you have generated all the figures and tables, otherwise the report may be incomplete.
```R
create_qc_reports("/path/to/sample/sheet", "screen", output_dir)
```

<p align="right">(<a href="#top">TOP</a>)</p>

<!-- Others -->
## Others
<a id="pandoc"></a>
### Pandoc
Pandoc is required to generate the R markdown report. Please download and install it from https://pandoc.org/installing.html

<p align="right">(<a href="#top">TOP</a>)</p>

<a id="test"></a>
### apeglm
the version of apeglm must be >= 1.22.1, optimHess problem in the lower version like below.
```
Error in optimHess(par = init, fn = nbinomFn, gr = nbinomGr, x = x, y = y,  : 
non-finite value supplied by optim
```

<p align="right">(<a href="#top">TOP</a>)</p>

<a id="conda"></a>
### Conda
When installing DESeq2, it may have error (Rlog1) on Mac M1 chip. Try cmd below to fix it.

```
export PKG_CPPFLAGS="-DHAVE_WORKING_LOG1P"
```

<p align="right">(<a href="#top">TOP</a>)</p>