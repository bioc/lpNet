\name{summarizeRepl}
\alias{summarizeRepl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarize Replicate Measurements
}
\description{
The function returns the the summarized replicate measuremenst.
}
\usage{
summarizeRepl(data,type=median)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The data matrix.
}
  \item{type}{
The summarization type which shall be used. Default: median.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Numeric matrix: the summarized data.
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
data("SahinRNAi2008")
## process data
dataStim <- dat.normalized[dat.normalized[,17]==1,-17]

# summarize replicates 
dataSt <- t(summarizeRepl(dataStim,type=mean))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{matrix summarization }


