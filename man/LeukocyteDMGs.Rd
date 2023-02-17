\name{LeukocyteDMGs}
\alias{LeukocyteDMGs}
\title{Identification of leukocyte differential genes in methylation beta value data from blood samples.}
\usage{
LeukocyteDMGs(platform,disease,normal,degene,freq,single_value,DMGs_value,FDR=TRUE)
}
\arguments{
   \item{platform}{Methylation data platform used(450 or 850).}
   \item{disease}{The methylation beta value dataframe of disease blood sample.The first columns is the gene symbols and the remaining columns are the gene methylation values of the disease samples.The gene symbols must be consistent with the normal.}
   \item{normal}{The methylation beta value dataframe of normal blood sample.The first columns is the gene symbols and the remaining columns are the gene methylation values of the normal samples.The gene symbols must be consistent with the disease.}
   \item{degene}{The differential methylation genes of normal and disease was identified using the T-test.}
   \item{frep}{The criteria for identifying stable gene pairs. The default setting of freq is 0.95.}
   \item{single_value}{The threshold for identifying individual-level differentially methylation genes.}
   \item{DMGs_value}{The threshold for identifying leukocyte differential methylation genes.}
   \item{FDR}{Whether to perform BH correction.The default setting is TRUE.}
}
\value{
A list with the following two items.

individual
A list with individual-level result.It contains all gene symbols,individual-level DMGs in each case sample (value is 1 for hypermethylated,value is -1 for hypomethylated) and adjusted pvalue
using the Benjamini and Hochberg approach.

population
A list with population-level DMGs.
}
\description{
LeukocyteDMGs could identify of leukocyte differential genes in methylation beta value data from blood samples.
}
\examples{
data(example)
DMGs=LeukocyteDMGs(450,disease,normal,degene,0.95,0.05,0.05,FDR = TRUE)
}
