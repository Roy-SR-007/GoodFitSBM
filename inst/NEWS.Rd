\name{NEWS for GoodFitSBM}
\title{NEWS for GoodFitSBM}
\encoding{UTF-8}

\section{Version 0.0.1}{\itemize{
\item This version of the \emph{GoodFitSBM} package deals with \emph{goodness-of-fit} test for a \emph{beta-SBM} (one of the three variants of SBMs discussed in \emph{Karwa et al. (2023)} [\url{https://doi.org/10.1093/jrsssb/qkad084}] used for modelling network data.
\item \code{goftest} is the primary function, which performs the Monte Carlo goodness-of-fit test for a particular Stochastic Blockmodel (SBM) viz., the beta SBM.
\item \code{sample_a_move} performs sampling of a graph; sampling done as per the beta SBM framework after one Markov move (basis).
\item \code{get_mle()} estimates the edge probabilities between blocks in a graph; using the maximum likelihood method of estimation, where the estimates are later required to compute the value of the chi-square test statistic (for performing the goodness-of-fit test) in \code{graphchi}.
\item \code{zachary} is the classic \emph{Zachary's Karate club data}, on which the goodness-of-fit test under a beta SBM framework is implemented.
}}
