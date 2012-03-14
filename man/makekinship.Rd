\name{makekinship}
\alias{makekinship}
\title{
  Create a kinship matrix
}
\description{
Compute the overall kinship matrix for a collection of families. This function and R documentation are included from the archived \code{kinship} package, which is 
available at http://cran.r-project.org/web/packages/kinship/index.html. 
}
\usage{
makekinship(famid, id, father.id, mother.id, unrelated=0)
}
\arguments{
\item{famid}{
a vector of family identifiers
}
\item{id}{
a vector of unique subject identifiers
}
\item{father.id }{
for each subject, the identifier of their biolgical father
}
\item{mother.id }{
for each subject, the identifier of thier biological mother
}
\item{unrelated}{
subjects with this family id are considered to be unrelated singletons,
i.e., not related to each other or to anyone else.
}
}
\value{
a kinship matrix
}
\examples{
\dontrun{
> ped <- read.csv("ped.csv")
> names(ped)
[1] "id"       "famid"    "fa"       "mo"       "sex"     
> kmat<-makekinship(ped$famid,ped$id,ped$fa,ped$mo)
}
}
