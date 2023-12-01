#' Zachary Karate Club Data
#'
#' This is the well-known and much-used Zachary karate club network. The data was collected from the members of a university karate club by Wayne Zachary in 1977. Each node represents a member of the club, and each edge represents a tie between two members of the club. The network is undirected. An often discussed problem using this dataset is to find the two groups of people into which the karate club split after an argument between two teachers.
#'
#' @source (Zachary, 1977),
#'   <http://vlado.fmf.uni-lj.si/pub/networks/data/Ucinet/UciData.htm#zachary>.
#' @format Two 34 by 34 matrices:
#' \describe{
#' \item{ZACHE}{symmetric, binary 34 by 34 adjacency matrix.}
#' \item{ZACHC}{symmetric, valued 34 by 34 matrix, indicating the relative strength of the associations}
#' }
"zachary"
