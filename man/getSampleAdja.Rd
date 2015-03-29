\name{getSampleAdja}
\alias{getSampleAdjaMAD}
\alias{getSampleAdja}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get The Sample Adjacency.
}
\description{
The function computes the adjacency of the edges computed in each step of the "loocv" or the "kfoldCV" function. If the variance of each edge shall be taken into account use "getSampleAdjaMAD", otherwise "getSampleAdja".
}
\usage{
getSampleAdjaMAD(edges_all, n, annot_node, method = median,
                 method2 = mad, septype = "->")
getSampleAdja(edges_all, n, annot_node, method = median, septype = "->")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{edges_all}{
The inferred edges using the "loocv" or the "kfoldCV" function.
}
  \item{n}{
Integer: the number of nodes.
}
  \item{annot_node}{
Vector of character strings: the annoation of the nodes.
}
  \item{method}{
Character string: the method used to summarize the edges of the individual steps. Default: "median".
}
  \item{method2}{
Character string: the method used for the computation of the variation of the edges of the individual steps. Default: "mad".
}
  \item{septype}{
Character string: the type of separation of two nodes in the annot string vector. Default: "->".
}
}

\value{
Numeric matrix: the adjacency matrix.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{loocv}}, \code{\link{kfoldCV}}
}
\examples{
# compute random edge weights
edges_all <- matrix(rnorm(5*6), nrow=5, ncol=6)

# annotation of the edges as returned by "loocv" and kfoldCV
colnames(edges_all) <- c("1->2", "1->3", "2->1", "2->3", "3->1", "3->2")

# annotation of the nodes
annot_node <- c(1,2,3)
getSampleAdjaMAD(edges_all, n=3, annot_node, method = "median", method2 = "mad", septype = "->")
  
getSampleAdja(edges_all, n=3, annot_node, method = "median", septype = "->")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{adjacency }

