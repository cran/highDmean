#' An example of GO term data
#'
#' A dataset containing the gene expressions for a Gene Ontology (GO) term
#' on two phenotype groups: BCR/ABL and NEG.
#' The id of the GO term is \code{GO:0000003}.
#' The raw dataset is taken from \code{ALL} package.
#' The data were preprocessed, for which the details are  elaborated in Zhang and Wang (2020).
#'
#'
#' @format A list with two subsets of gene expression data.
#' \describe{
#'   \item{X}{A matrix containing gene expressions for the BCR/ABL group. The row id is for patient and the
#'   column id is for gene.}
#'   \item{Y}{A matrix containing gene expressions for the NEG group. The row id is for patient and the
#'   column id is for gene.}
#' }
#' @references
#' Zhang, H. and Wang, H. (2020). Result consistency of high dimensional
#' two-sample tests applied to gene ontology terms with gene sets. Manuscript in review.
#'
"GO_example"
