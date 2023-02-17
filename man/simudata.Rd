\name{simudata}
\alias{simudata}
\title{Generate data for simulation experiments.}
\usage{
simudata(platform,mye_ratio,lym_ratio)
}
\arguments{
   \item{platform}{Select the simulated methylation platform(450 or 850).}
   \item{mye_ratio}{Proportion of myeloid cells in simulated data.}
   \item{lym_ratio}{Proportion of lymphocytes in simulated data.}
}
\value{
A list with the following three items.

DMG
Positive DMG gene set in simulated data.

simu.normal
A beta value of dataframe for simulated normal blood.

simu.disease
A beta value of dataframe for simulated disease blood.
}
\description{
The code of simudata could generate data for simulation experiments.
}
\examples{
simu_0.6 = simudata(450,0.6,0.4)
}
