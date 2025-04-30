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
8. [For Developers](#for-developers) 

</details>

<!-- Dependencies-->
## Dependencies
```R
install.packages("configr", version = '0.3.5')
install.packages("vroom", version = '1.6.3')
install.packages("data.table", version = '1.14.8')
install.packages("Ckmeans.1d.dp", version = '4.3.5')
install.packages("gplots", version = '3.1.3')
install.packages("ggplot2", version = '3.4.3')
install.packages("plotly", version = '4.10.4')
install.packages("ggcorrplot", version = '0.1.4.1')
install.packages("corrplot", version = '0.92')
install.packages("see", version = '0.8.0')
install.packages("ggbeeswarm", version = '0.7.2')
install.packages("reactable", version = '0.4.4')
install.packages("reshape2", version = '1.4.4')
install.packages("htmltools", version = '0.5.6')
install.packages("sparkline", version = '2.0')
install.packages("dendextend", version = '1.17.1')
install.packages("gtools", version = '3.9.4')

install.packages("BiocManager")
BiocManager::install(version = "3.18")
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
library(sparkline)
library(dendextend)
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
| sample_name | replicate | condition | ref_time_point | library_independent_count | library_dependent_count | valiant_meta | vep_anno | adapt5 | adapt3 | per_r1_adaptor | per_r2_adaptor | library_name | library_type|
| - | - | - | - | - | - | - | - | - | - | - | - | - | - |
| sample1 | R1 | D4 | D4 | s1.allcounts.tsv.gz | s1.libcounts.tsv.gz | meta.csv.gz | meta_consequences.tsv.gz | CTGACTGGCACCTCTTCCCCCAGGA | CCCCGACCCCTCCCCAGCGTGAATG | 0.21 | 0.10 | libA | screen |
| sample2 | R2 | D4 | D4 |s2.allcounts.tsv.gz | s2.libcounts.tsv.gz | meta.csv.gz | meta_consequences.tsv.gz | CTGACTGGCACCTCTTCCCCCAGGA | CCCCGACCCCTCCCCAGCGTGAATG | 0.11 | 0.02 | libA | screen |
| sample3 | R3 | D4 | D4 |s3.allcounts.tsv.gz | s3.libcounts.tsv.gz | meta.csv.gz | meta_consequences.tsv.gz | CTGACTGGCACCTCTTCCCCCAGGA | CCCCGACCCCTCCCCAGCGTGAATG | 0.01 | 0.18 | libA | screen |
| sample4 | R1 | D7 | D4 | s4.allcounts.tsv.gz | s4.libcounts.tsv.gz | meta.csv.gz | meta_consequences.tsv.gz | CTGACTGGCACCTCTTCCCCCAGGA | CCCCGACCCCTCCCCAGCGTGAATG | 0.21 | 0.10 | libA | screen |
| sample5 | R2 | D7 | D4 |s5.allcounts.tsv.gz | s5.libcounts.tsv.gz | meta.csv.gz | meta_consequences.tsv.gz | CTGACTGGCACCTCTTCCCCCAGGA | CCCCGACCCCTCCCCAGCGTGAATG | 0.11 | 0.02 | libA | screen |
| sample6 | R3 | D7 | D4 |s6.allcounts.tsv.gz | s6.libcounts.tsv.gz | meta.csv.gz | meta_consequences.tsv.gz | CTGACTGGCACCTCTTCCCCCAGGA | CCCCGACCCCTCCCCAGCGTGAATG | 0.01 | 0.18 | libA | screen |

* *please use the same headers in the example*
* *replicate, condition and ref_time_point are optional, but required for screen qc*
* *adapt5 and adapt3 are optional, please leave them blank if you don't have them, but required for reads with primers*
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
Reference samples must be assigned. MAVEQC automatically creates reference samples (```maveqc_ref_time_point``` and ```maveqc_ref_time_point_samples```) from the sample sheet using ```ref_time_point``` and ```sampe_name```.

```R
output_dir <- "/path/to/output/directory"

samqc <- create_sampleqc_object(sge_objs)
samqc <- run_sample_qc(samqc, "screen")

qcplot_samqc_all(samqc, qc_type = "screen", plot_dir = output_dir)
qcout_samqc_all(samqc, qc_type = "screen", out_dir = output_dir)
```

<p align="right">(<a href="#top">TOP</a>)</p>

<a id="sqc2"></a>
### QC 2: Experimental QC
MAVEQC automatically creates the coldata (```maveqc_deseq_coldata```) from sample sheet for Screen QC.

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
expqc <- create_experimentqc_object(samqc) 
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


## For Developers

### Development Environment Setup

To quickly set up a development environment for this project, you can use Docker Compose. 
This will create all the necessary containers and configurations for a smooth development experience.

#### Requirements

Before running `docker compose up`, make sure you have the following installed: Docker and Docker Compose.

#### Running the Development Environment

Once you have Docker and Docker Compose installed, follow these steps:

1. **Clone the repository**:
```
git clone https://gitlab.internal.sanger.ac.uk/sci/MAVEQC.git
cd MAVEQC
```

2. **Build and start the development containers:**

```
docker compose up
```
This command will build the Docker image and start the MAVEQC development container (maveqc-dev-container) defined in 
docker-compose.yml. This will also mount the MAVEQC folder onto the container.

If you want to run the container in detached mode, use the -d flag:
```
docker compose up -d
```
This starts the container in the background, freeing the terminal for other commands.

3. **Access the maveqc-dev-container**
```
docker exec -it maveqc-dev-container bash
```
4. **Run MAVEQC**

Once inside the maveqc-dev-container, start R and load all the package components:
```
R
```
Then, load all the package components:
```
devtools::load_all()
```
After running devtools::load_all(), refer to the documentation above for instructions on using the MAVEQC reporter.
> **Note**: Use paths relative to `/usr/src/MAVEQC-R`, the container’s working directory.

To stop the running container, use:

```
docker compose down
```

#### Pull docker image from GitLab Container Registry
You can pull latest image from `gitlab-registry.internal.sanger.ac.uk` using the following steps:

1. Login into GitLab container registry using username and password:

```
docker login gitlab-registry.internal.sanger.ac.uk
```

2. Then, pull docker image from container registry using the following command:

```
docker pull gitlab-registry.internal.sanger.ac.uk/sci/maveqc
```

3. Once you pull the docker images make sure you replace exisiting `image` name in `docker-compose.yml` to `gitlab-registry.internal.sanger.ac.uk/sci/maveqc:latest`, you can see the following example:

```
services:
  maveqc-dev:
    image: gitlab-registry.internal.sanger.ac.uk/sci/maveqc:latest
...
```

4. Finally, you run docker container using docker compose command using the following command:

```
docker compose up -d
``` 
