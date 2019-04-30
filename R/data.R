# MS-EmpiRe - Mass Spectrometry analysis using Empirical and Replicate based statistics
# Copyright (C) 2018  Markus Gruber
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


#' @import Biobase

#' @export
read.standard <- function(f, sample.mapping, signal_pattern="c", sep="\t", id_col = 1, prot.id.col=NULL,
                          prot.id.generator=NULL, remove.pattern=F, comment.char="#")
{
  if((!is.null(prot.id.col)) && (!is.null(prot.id.generator)))
  {
    stop("both prot.id.generator and prot.id.col were not NULL. Please enter only one of theses two arguments")
  }
  
  message(sprintf("Reading data from %s", f))
  
  data <- read.csv(f, sep=sep, stringsAsFactors = FALSE, comment.char=comment.char)

  data[, id_col] <- as.character(data[, id_col])
  
  # reorder data by prot.id
  if(!is.null(prot.id.generator))
  {
    prot.ids <- sapply(data[, id_col], prot.id.generator)
    data <- data[order(prot.ids), ]
  }else if(!is.null(prot.id.col)){
    prot.ids <- data[, prot.id.col]
    data <- data[order(prot.ids), ]
  }

  ids <- data[,id_col]

  signal_cols <- grep(signal_pattern,colnames(data))

  exprs <- as.matrix(apply(data[, signal_cols], 2, as.numeric))

  if(remove.pattern==T)
  {
    colnames(exprs) <- gsub(signal_pattern, "", colnames(exprs))
  }
  rownames(exprs) <- ids


  f_data <- as.data.frame(data[, setdiff(1:ncol(data), c(signal_cols, id_col))], stringsAsFactors=FALSE)
  colnames(f_data) <- colnames(data)[setdiff(1:ncol(data), c(signal_cols, id_col))]
  rownames(f_data) <- ids


  p_data <- read.csv(sample.mapping, sep=sep, header=T, row.names = 1, stringsAsFactors = FALSE)
  p_data <- as.data.frame(p_data[colnames(exprs), ], stringsAsFactors=FALSE)
  rownames(p_data) <- colnames(exprs)
  colnames(p_data) <- c("condition")
  p_data <- Biobase::AnnotatedDataFrame(p_data)

  res <- Biobase::ExpressionSet(exprs, p_data)

  if(!is.null(prot.id.col))
  {
    f_data["prot.id"] <- as.character(f_data[, prot.id.col])
  }
  if(!is.null(prot.id.generator))
  {
    prot_col <- sapply(rownames(exprs), prot.id.generator)
    f_data["prot.id"] <- as.character(prot_col)
  }

  fData(res) <- f_data
  return(res)
}

#' @export
read.EB <- function(dir, expression="exprs.txt", pheno="p_data.txt", feature="f_data.txt", sep="\t", header=F)
{
  
  exprs <- as.matrix(read.csv(file.path(dir, expression), sep="\t", header=header, stringsAsFactors = FALSE))
  
  pdat <- read.csv(file.path(dir, pheno), sep="\t", stringsAsFactors = FALSE, header=header)
  if(ncol(pdat) != 2)
  {
    stop("pheno data file should only contain two columns, sample and condition!")
  }
  if(nrow(pdat) != ncol(exprs))
  {
    stop("number of samples in pheno data file and exprs file are not the same!")
  }
  
  fdat <- read.csv(file.path(dir, feature), sep="\t", header=header, stringsAsFactors = FALSE)
  
  if(nrow(fdat) != nrow(exprs))
  {
    stop("number of features in feature data file and exprs file are not the same!")
  }
  
  colnames(pdat) <- c("sample", "condition")
  rownames(pdat) <- pdat$sample
  
  rownames(fdat) <- fdat[, 1]
  fdat$prot.id <- rownames(fdat)
  
  colnames(exprs) <- rownames(pdat)
  rownames(exprs) <- rownames(fdat)
  
  res <- Biobase::ExpressionSet(exprs, Biobase::AnnotatedDataFrame(pdat))
  fData(res) <- fdat
  return(res)
}




#' @export
read.MaxQuant <- function(peptides, sample.mapping, sep="\t", feature.names="id")
{
  return(read.standard(peptides,
                       sample.mapping=sample.mapping,
                       signal_pattern="Intensity\\.",
                       remove.pattern=T,
                       id_col=feature.names
                       ))
}
