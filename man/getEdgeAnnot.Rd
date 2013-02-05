\name{getEdgeAnnot}
\alias{getEdgeAnnot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get the annotation of the edges.
}
\description{
The function returns the annotation of the edges needed for the LP. Positive edges are annotated with "w+" and negative with "w-". The given nodes are just enumerated from 1 to n and the edge between node i and j is given by "w+_i_j" for the positive, respectively by "w-_i_j" for the negative edges. The annotation "w_i_^_0" defines the baseline activity of gene i.
}
\usage{
getEdgeAnnot(n,allpos)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Integer: number of genes.
}
  \item{allpos}{
Logical: should all edges be positive? Corresponds to learning only activating edges. Default: FALSE.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Vector of character strings.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{

}
\examples{
n <- 5
annot <- getEdgeAnnot(n)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{annotation }


