# LeukocyteDMGs
LeukocyteDMGs is an algorithm to identify of leukocyte differential genes in methylation array data.
# Install
To install the LeukocyteDMGs, install from github using devtools:
```
library(devtools)
install_github("FunMoy/LeukocyteDMGs")
```
# Usage
LeukocyteDMGs(platform,disease,normal,degene,freq,single_fdr,DMGs_fdr)
Arguments|Description
:--|:---
platform|Methylation data platform used(450 or 850).
disease|The dataframe of disease sample.The first columns is the gene symbols and the remaining columns are the gene methylation values of the disease samples.The gene symbols must be consistent with the normal.
normal|The dataframe of normal sample.The first columns is the gene symbols and the remaining columns are the gene methylation values of the normal samples.The gene symbols must be consistent with the disease.
degene|The differential methylation genes of normal and disease was identified using the T-test.
frep|The criteria for identifying stable gene pairs. The default setting of freq is 0.95.
single_fdr|The threshold of FDR for identifying individual-level differentially methylation genes.
DMGs_fdr|The threshold of FDR for identifying leukocyte differential methylation genes.}



# Example
```
library(LeukocyteDMGs)
data(example)
DMGs=LeukocyteDMGs(450,disease,normal,degene,0.95,0.05,0.05)
```

# Contact email
Please don't hesitate to address comments/questions/suggestions regarding this R package to:
Qi Fan <FunMoy@163.com>; Haidan Yan <Joyan168@126.com>
