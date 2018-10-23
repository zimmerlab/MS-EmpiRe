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
read.standard <- function(f, sample.mapping=NULL, signal_pattern="c", sep="\t", id_col = 1, prot.id.col=NULL,
                          prot.id.generator=NULL, remove.pattern=F)
{
  message(sprintf("Reading data from %s", f))
  data <- read.csv(f, sep=sep, stringsAsFactors = FALSE)

  
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


  if(remove.pattern==T)
  {
    colnames(exprs) <- gsub(signal_pattern, "", colnames(exprs))
  }
  rownames(exprs) <- ids


  f_data <- as.data.frame(data[, setdiff(1:ncol(data), c(signal_cols, id_col))])
  colnames(f_data) <- colnames(data)[setdiff(1:ncol(data), c(signal_cols, id_col))]
  rownames(f_data) <- ids


  p_data <- read.csv(sample.mapping, sep=sep, header=T, row.names = 1, stringsAsFactors = FALSE)
  p_data <- as.data.frame(p_data[colnames(exprs), ])
  rownames(p_data) <- colnames(exprs)
  colnames(p_data) <- c("condition")
  p_data <- Biobase::AnnotatedDataFrame(p_data)

  res <- Biobase::ExpressionSet(exprs, p_data)

  if(!is.null(prot.id.col))
  {
    res["prot.id"] <- as.character(res[prot.id.col])
    res[prot.id.col] <- NULL
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
read.EB <- function(dir, expression="exprs.txt", pheno="p_data.txt", feature="f_data.txt", sep="\t")
{
  data <- Biobase::readExpressionSet(exprsFile = file.path(dir, expression),
                     phenoDataFile = file.path(dir, pheno),
                     sep=sep)
  fData(data) <- read.csv(file.path(dir, feature), sep="\t", stringsAsFactors = FALSE)

  tmp <- fData(data)
  rownames(tmp) <-  tmp[,1]
  tmp[,1] <- NULL
  fData(data) <- tmp
  return(data)
}

#' @export
read.MaxQuant <- function(peptides, sep="\t", sample.mapping=NULL, feature.names="id")
{
  return(read.standard(peptides,
                       sample.mapping=sample.mapping,
                       signal_pattern="Intensity\\.",
                       remove.pattern=T,
                       id_col=feature.names
                       ))
}
