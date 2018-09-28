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

#' @param data MSnSet
#' @param rate minimal amount of samples per condition in which each peptide must be detectable
#' @return MSnSet with the filtered exprs and fData fields
#' @export
filter_detection_rate <- function(data, rate=2, condition=NULL)
{
  x <- exprs(data)

  if(is.null(condition))
  {
    condition <- extract_conditions(data)
  }

  s <- (x > 0) %*% condition

  f <- apply(s >= rate, 1, all)

  print(sprintf("Removed %d peptides because of low detection rate, i.e. detected in less than %d samples)", nrow(x) - sum(f), rate))

  samples_to_keep <- rowSums(condition) > 0

  exp <- exprs(data)[f, samples_to_keep]
  f_dat <- fData(data)[f, ]
  p_dat <- data.frame(condition=pData(data)$condition[samples_to_keep], row.names=rownames(pData(data))[samples_to_keep])
  p_dat <- Biobase::AnnotatedDataFrame(p_dat)
  r <- Biobase::ExpressionSet(exp, p_dat)
  fData(r) <- f_dat
  return(r)
}

check_sign <- function(s, sign)
{
  return(s == sign)
}

#' @export
filter_MaxQuant <- function(data, proteinGroups, filter_on=c("Potential.contaminant", "Reverse"), filter_symbol="+",
                            extern.filter=NULL, valid.peptides=NULL)
{
  proteinGroups <- read.csv(proteinGroups, sep="\t")


  to_filter <- sapply(filter_on, function(f)
    {
    check_sign(proteinGroups[, f], filter_symbol)
  })

  #contains T/F whether protein should be kept for analysis
  to_filter <- apply(to_filter, 1, any)

  if(!is.null(extern.filter))
  {
    extern <- sapply(proteinGroups$Protein.IDs, extern.filter)
    print(sprintf("extern: %d", sum(extern)))
    to_filter <- to_filter | extern
  }


  print(sprintf("will remove %d proteins from analysis", sum(to_filter)))


  good_prots <- proteinGroups[!to_filter, c("id", "Protein.IDs", "Peptide.IDs")]
  print(sprintf("remaining for analysis: %d", nrow(good_prots)))


  good_prots <- apply(good_prots, 1, function(prot)
    {
    peps <- unlist(strsplit(prot[3], ";"))
    m <- matrix(c(peps, rep(prot[2], length(peps))), ncol=2, nrow=length(peps))
    return(m)
  })

  good_prots <- do.call("rbind", good_prots)
  colnames(good_prots) <- c("pep.id", "prot.id")

  tmp_exprs <- exprs(data)[good_prots[, "pep.id"], ]
  tmp_f_data <- fData(data)[good_prots[,"pep.id"], ]
  tmp_p_data <- pData(data)


  rnames <- apply(good_prots, 1, function(r) paste(c(r[2], r[1]), collapse="."))
  rownames(tmp_f_data) <- rnames
  rownames(tmp_exprs) <- rnames

  tmp_f_data$prot.id <- good_prots[, "prot.id"]

  if(!is.null(valid.peptides))
  {
    keep <- tmp_f_data[, "id"] %in% valid.peptides
    tmp_f_data <- tmp_f_data[keep, ]
    tmp_exprs <- tmp_exprs[keep, ]
  }

  r <- Biobase::ExpressionSet(tmp_exprs, Biobase::AnnotatedDataFrame(tmp_p_data))
  fData(r) <- tmp_f_data
  return(r)
}
