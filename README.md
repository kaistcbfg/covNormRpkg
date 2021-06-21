# covNormRpkg 

## Introduction

Hi-C and capture Hi-C have greatly advanced our understanding the principle of higher order chromatin structure, but still development of computational methods that effectively eliminate innate biases is essential to precise detection of chromatin contacts. To resolve this issue, we developed an implicit normalization method named ‘covNorm’ and implemented as an R package which performs negative binomial model-based normalization for genomic bin coverages and ligated DNA fragment distances, visualizes the process for quality control, and provides significance scores of the chromatin contacts. Taken together, the proposed method provides accurate and reproducible results between biological replicates as well as easily applicable to various ‘C’ technologies.

## Publication and Citation
If you use covNorm in a paper, please cite:

> Kyukwang Kim and Inkyung Jung,  
> covNorm: an R package for coverage based  normalization of Hi-C and capture Hi-C data,  
> *Computational and Structural Biotechnology Journal*, Volume 19, pages 3149-3159(2021).  
> doi: https://doi.org/10.1016/j.csbj.2021.05.041 

The preliminary forms of covNorm in previous works:

> Kim, K., *et al*,
> 3DIV update for 2021: a comprehensive resource of 3D genome and 3D cancer genome  
> *Nucleic Acids Research*, Volume 49, Issue D1, pages D38–D46(2020).  
> doi: https://doi.org/10.1093/nar/gkaa1078

> Jung, I., *et al*,   
> A compendium of promoter-centered long-range chromatin interactions in the human genome   
> *Nature Genetics* Volume 51, pages 1442–1449(2019).  
> doi: https://doi.org/10.1038/s41588-019-0494-8

> Yang, D., *et al*,   
> 3DIV: A 3D-genome Interaction Viewer and database   
> *Nucleic Acids Research*, Volume 46, Issue D1, pages D52–D57(2018).  
> doi: https://doi.org/10.1093/nar/gkx1017

## License
Copyright (c) YEAR: 2020 COPYRIGHT HOLDER: KAIST (Corporation registration number: 114471-0000668).
Registered at the Korea Copyright Commission(C-2021-022800) in accordance with Article 53 of the Copyright Act. 

Developed by Kyukwang Kim & Inkyung Jung, KAIST Dept. of Biological Sciences.
For commercial use of the software, please contact the authors.

## Installation

Prerequisites for covNormRpkg can be installed by
```R
install.packages(c("MASS", "propagate", "FAdist", "stringr", "splines"))    # Imports
install.packages(c("reshape2", "gplots", "ggplot2", "corrplot"))            # Suggests
```
Install the latest released version from GitHub using [devtools](https://cran.r-project.org/package=devtools) with

```R
devtools::install_github("kaistcbfg/covNormRpkg")
```
Download covNorm source code and install with the R CMD

```bash
git clone https://github.com/kaistcbfg/covNormRpkg.git
R CMD build covNormRpkg
R CMD INSTALL covNormRpkg_1.1.0.tar.gz
```

## Sample Data and Input Format

Sample Hi-C/pcHi-C datasets for the covNorm can be downloaded with:

```bash
wget http://junglab.kaist.ac.kr/Dataset/GM19204.chr17.cis.feature.gz #Hi-C Replicate 1,  783588 rows, 12Mb
wget http://junglab.kaist.ac.kr/Dataset/GM19240.chr17.cis.feature.gz #Hi-C Replicate 2, 1136435 rows, 19Mb

wget http://junglab.kaist.ac.kr/Dataset/GM12878.po.no0.feature.gz #pcHi-C Replicate 1, 4298847 rows, 59Mb
wget http://junglab.kaist.ac.kr/Dataset/GM19240.po.no0.feature.gz #pcHi-C Replicate 2, 3582541 rows, 49Mb
```  

The input file must follow following format, imported as a data.frame in R.
|        frag1        |           frag2          | cov_frag1 | cov_frag2 | freq | dist     |
|:-------------------:|:------------------------:|:---------:|:---------:|------|----------|
| chr17.140000.160000 |  chr17.83160000.83180000 |    2296   |    2304   |  1.0 | 83020000 |
| chr17.140000.160000 | chr17.83180000.83200000  |    2296   |    2072   |  2.0 | 83040000 |
| chr17.140000.160000 |  chr17.83200000.83220000 |    2296   |    778    |  2.0 | 83060000 |
|         ...         |            ...           |    ...    |    ...    |  ... |    ...   |
| chr17.160000.180000 |    chr17.200000.220000   |    2119   |    2253   | 12.0 |   40000  |
| chr17.160000.180000 |    chr17.220000.240000   |    2119   |    1744   |  9.0 |   60000  |

**frag1**: Dot spliced chromosome, start coordinate, end coordinate of the first DNA fragment.  
**frag2**: Dot spliced chromosome, start coordinate, end coordinate of the second DNA fragment.  
In case of pcHi-C, frag1 should be 'promoter' site in PO-interaction normalization. 

**cov_frag1**: Coverage value of frag1 bin.  
**cov_frag2**: Coverage value of frag2 bin.  
Use 'coverageBed' tool from ['bedtools'](https://github.com/arq5x/bedtools2) to compute bin coverage.

**freq**: Raw interaction frequency between two bins.   
**dist**: Genomic distance between two bins' start coordinates.  
Use midpoint between start and end for pcHi-C.


## Example Codes

Please refer to the 'Supplmentary Information' of our published paper or the man page of each function (e.g. ```?covNormRpkg::normCoverage```) for analysis tips and details about input parameters.

We  provide example codes for Hi-C and pcHi-C analysis demonstration. 

```R
library('covNormRpkg')

args <- commandArgs(TRUE)
file_name=args[1]

raw_data <- read.table(gzfile(file_name),head=TRUE)
print("1: Data Loaded.")

raw_data_filter <- covNormRpkg::filterInputDF(raw_data)
print("2: Data Filtered.")

cov_result <- covNormRpkg::normCoverage(raw_data_filter)
cov_result$coeff_cov1
cov_result$coeff_cov2
cov_df <- cov_result$result_df
write.table(cov_df, file=gzfile("outFileName1"), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE) # Coverage normalization result
print("3: Coverage normalized.")

covNormRpkg::checkFreqCovPCC(cov_df, outpdfname='QCplot_coverage_PCC.pdf')
covNormRpkg::plotCovNormRes( cov_df, outpdfname='QCplot_coverage_heatmap.pdf')
print("4: Plot coverage normalization results.")

dist_result <- covNormRpkg::normDistance(cov_df, max_dist=2000000)
dist_df <- dist_result$result_df
print("5: Distance normalized.")

covNormRpkg::checkFreqDistPCC(dist_df, outpdfname='QCplot_dist_PCC.pdf')
covNormRpkg::plotDistNormRes( dist_df, outpdfname='QCplot_dist_hexmap.pdf')
print("6: Plot distance normalization results.")

final_df <- covNormRpkg::contactPval(dist_df, 'fit.pdf')
print("7: Significant interactions called.")

#Uncomment 'saveEachChr' to split-save file for each chromosome.
#covNormRpkg::saveEachChr(final_df, "./outputFolder", "outputSampleName") 
write.table(final_df, file=gzfile("outFileName2"), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE) # Distance normalization & significant interactions
```
For pcHi-C PO-interaction normalization, modify step 3 (Coverage normalization) part with following code.
```R
raw_data_filter <- raw_data_filter[which(raw_data_filter$dist<2000000),] # Consider 15kb-2Mb interactions only in pcHi-C
cov_result <- covNormRpkg::normCoverage(raw_data_filter, do_shuffle=FALSE, cov1_thresh=200, cov2_thresh=50) 
#in pcHi-C PO-interaction normalization, lower coverage threshold for 'other' interaction
```

## Output 

For detailed interpretation of the result, check 'Supplmentary Information' of our published paper.

Following columns will be added to the input dataframe:

| rand | exp_value_capture | capture_res | exp_value_dist | dist_res | p_result_dist | FDR_dist_res |
|:----:|:-----------------:|:-----------:|:--------------:|:--------:|:-------------:|:------------:|
|  77  |       2.0625      |    2.9091   |     1.7021     |  1.4464  |    0.053347   |   0.709688   |
|  47  |       1.049       |    1.9066   |     1.6201     |  1.1093  |    0.237543   |   0.945149   |
|  40  |       1.6197      |    0.6174   |     0.7307     |  0.9346  |    0.475677   |   0.948673   |
|  38  |       1.7874      |    0.5595   |     0.6794     |  0.9287  |    0.486306   |   0.948673   |
|  ... |        ...        |     ...     |       ...      |    ...   |      ...      |      ...     |
|  84  |       0.7842      |    1.2752   |     1.3979     |  0.9489  |    0.450677   |   0.948673   |

**rand**: Random integer used for coverage 1 & 2 shuffle. Not used in pcHi-C normalization.  
**exp_value_capture**: Expected interaction frequency between DNA fragment based on coverage values.  
**capture_res**: Coverage normalized interaction frequency between DNA fragments.  
**exp_value_dist**: Expected interaction frequency between DNA fragment based on genomic distance between DNA fragments.  
**dist_res**: Distance normalized interaction frequency between DNA fragments.   
**p_result_dist**: *p*-value of the interaction.  
**FDR_dist_res**: FDR of the interaction.  


Use **'capture_res'** value to plot coverge normalized Hi-C contact map.  
Following example code can be used:

```python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import gzip

chromosome_size = 83257441 #hg38 chr17
resolution = 40000 # 40kb
bin_length = (chromosome_size/resolution) + 1

contact_map = np.zeros((bin_length, bin_length))

f = gzip.open('GM19240.chr17.covnorm.gz') # coverage normalization result file
f.readline() #header
for line in f:
    line = line.rstrip()
    linedata = line.split('\t')
    bin1 = int(linedata[0].split('.')[1])/resolution
    bin2 = int(linedata[1].split('.')[1])/resolution
    freq = float(linedata[8])

    contact_map[bin1][bin2] += freq
    contact_map[bin2][bin1] += freq
#
f.close()

fig = plt.figure(1)
ax = fig.add_subplot(111)
cax = ax.matshow(contact_map, cmap=plt.cm.RdYlBu_r, vmin=0, vmax=5)
fig.colorbar(cax)
plt.savefig("HiC_contact_map.pdf", dpi=1000)
```
GM19240 chr17 normalized Hi-C contact matrix  

<img src="https://user-images.githubusercontent.com/67453667/87240706-f4b5ad00-c456-11ea-9ae3-8263be2e2ebd.png" width="50%"></img>
