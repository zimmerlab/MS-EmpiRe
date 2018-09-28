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



generate_foldchanges <-function(x, cond1, cond2)
{
  c1 <- x[, cond1]
  c2 <- x[, cond2]
  pairs <- t(expand.grid(1:sum(cond2), 1:sum(cond1)))
  # pairs <- t(expand.grid(1:sum(cond1), 1:sum(cond2)))
  pairs <- pairs[c(2,1), ]

  res <- apply(pairs, 2, function(pair) log2(c2[, pair[2]] / c1[, pair[1]]))
  return(res)
}



norm_background <- function(diffs, diffs2bg, signs, background_distributions, background_zeros)
{
  return(sapply(1:length(diffs), function(pep)
    {
    zscore <- diffs[pep]
    if(zscore <= 0)
    {
      return(0.0)
    }
    bg_dist <- background_distributions[[diffs2bg[pep]]]

    where <- (length(bg_dist) + 1) - bin_search(bg_dist, zscore) + 1

    rank <- where / (length(bg_dist) + 1)

    rank <- qnorm(max(0.5, 1.0 - rank))
    rank <- rank * signs[pep]
    return(rank)
  }))
}


toZ_ranks <- function(peps, fcs, assignment, dists)
{
  return(apply(cbind(peps, 1:length(peps)), 1, function(r)
  {
    peptide <- r[1]
    peptide_loc_idx <- r[2]

    idx <- assignment[peptide]
    distrib <- dists[[idx]]

    loc_fcs <- fcs[peptide_loc_idx, ]
    ranks <- rep(-1, length(loc_fcs))
    finite <- is.finite(loc_fcs)


    loc_fcs <- loc_fcs[is.finite(loc_fcs)]
    positive <- loc_fcs > 0

    tmp_ranks <- sapply(loc_fcs, function(fc) bin_search(distrib, fc)) # returns smaller / eq including fc

    if(sum(positive) > 0)
    {
      tmp_ranks[positive] <- (length(distrib) + 1) - tmp_ranks[positive] + 1
    }

    ranks[finite] <- tmp_ranks
    return(ranks)
  }))
}

toZ <- function(peps, fcs, assignment, dists)
{
  return(apply(cbind(peps, 1:length(peps)), 1, function(r)
    {
      peptide <- r[1]
      peptide_loc_idx <- r[2]

      MAX_Z = 7.650628
      MAX_Z_P = 1.0 - 10e-15

      idx <- assignment[peptide]
      distrib <- dists[[idx]]

      loc_fcs <- fcs[peptide_loc_idx, ]

      loc_fcs <- loc_fcs[is.finite(loc_fcs)]
      positive <- loc_fcs > 0
      positive2 <- loc_fcs >= 0


      ranks <- sapply(loc_fcs, function(fc) bin_search(distrib, fc)) # returns smaller / eq including fc

      if(sum(positive) > 0)
      {
        ranks[positive] <- (length(distrib) + 1) - ranks[positive] + 1
      }

      ranks <- ranks / (length(distrib) + 1)
      ranks <- abs(qnorm(ranks))
      if(sum(!positive2) > 0)
      {
        ranks[!positive2] <- -ranks[!positive2]
      }

      # each of the ranks is in N(0,1) space
      return(sum(ranks))
  }))
}

#' @export
de.ana <- function(data, out.dir=NULL, with_fcs=F)
{
  message("starting differential analysis...")
  #build condition2sample mapping
  conditions <- extract_conditions(data)

  cond1 <- conditions[,1] == 1
  cond2 <- ! cond1

  #generate needed pairs for background (3v3, 1v2, 2v4...)
  filter <- exprs(data) != 0
  ttt <- cbind(rowSums(filter[, cond1]), rowSums(filter[, cond2]))
  rownames(ttt) <- NULL
  ttt <- t(apply(ttt, 1, sort))
  pep2validpair <- apply(ttt, 1, function(r) paste(r, collapse="v"))
  pep2fullvar <- apply(ttt, 1, function(r) num_covar(r) + num_var(r))


  uniq_rep_pairs <- unique(ttt)


  background_distributions <-lapply(split(uniq_rep_pairs, seq(nrow(uniq_rep_pairs))), function(w) generate_distribution(w, out.dir=out.dir))
  names(background_distributions) <- apply(uniq_rep_pairs, 1, function(r) paste(r, collapse="v"))
  background_zeros <- lapply(background_distributions, function(d) sum(d==0))


  #build contexts
  ctx <- msEmpiRe::build_context(data, cond1, cond2, out.dir=out.dir)
  dists <- ctx[[1]]
  assignment <- ctx[[2]]
  sdsss <- ctx[[3]]


  if(!is.null(out.dir))
  {
    pdf(file.path(out.dir, "bg_sds.pdf"))
    plot(sdsss)
    dev.off()
  }
  pep2sds <- sapply(1:nrow(exprs(data)), function(i) sdsss[assignment[i]])


  #actual analysis
  xprs <- exprs(data)
  fcs <- generate_foldchanges(xprs, cond1, cond2)

  tmp <- cbind(fcs, pep2sds)


  random_fcs <- t(apply(tmp, 1, function(r)
  {
    return(to_normal(r[1:ncol(fcs)], r[ncol(fcs)+1]))
  }))


  what <- aggregate(1:nrow(fData(data)) ~ prot.id, fData(data), c)
  num_peps <- aggregate(1:nrow(fData(data)) ~ prot.id, fData(data), length)
  colnames(what) <- c("prot.id", "pep.ids")

  FC_LIMITS <- c(0.5, 0.65, 0.8, 0.9, 0.95)

  message("computing peptide level scores...")

  test_scores <- cbind(toZ(1:nrow(fcs), fcs, assignment, dists), pep2fullvar)
  test_scores <- cbind(test_scores,
                       sapply(1:nrow(fcs), function(ii) sdsss[assignment[ii]]),
                       sapply(1:nrow(fcs), function(ii) assignment[ii]))

  medians <- apply(random_fcs, 1, function(r)
  {
    ttt <- r[is.finite(r)]
    ttt <- sort(ttt)
    return(ttt[0.5*length(ttt) + 1])
  }
    )
  test_scores <- cbind(test_scores, medians)

  fcs_shifted <- fcs - medians

  shifted_scores <- toZ(1:nrow(fcs_shifted), fcs_shifted, assignment, dists)

  test_scores <- cbind(test_scores, shifted_scores)

  diffs <- abs(test_scores[, 1]) - abs(shifted_scores)
  s <-sign(test_scores[, 1])
  s[s==0] <- 1

  ranks <- t(toZ_ranks(1:nrow(fcs), fcs_shifted, assignment, dists))
  test_scores <- cbind(test_scores, ranks)
  test_scores <- cbind(test_scores, fcs)

  bg_normed_score <- norm_background(diffs, pep2validpair, s, background_distributions, background_zeros)

  test_scores <- cbind(test_scores, bg_normed_score)
  test_scores <- cbind(test_scores, diffs)
  test_scores <- cbind(test_scores, fcs_shifted)
  test_scores <- cbind(test_scores, sapply(1:nrow(fcs), function(pep) length(dists[[assignment[pep]]]))  )
  rownames(test_scores) <- rownames(exprs(data))
  colnames(test_scores) <- c("score", "sd", "bg_sd", "bg_idx", "median", "bg_score", paste0("rank", 1:ncol(ranks)), paste0("FC", 1:ncol(fcs)), "bg_normed_score", "difference", paste0("FC_SH", 1:ncol(fcs)), "ctx_size")

  message("merging scores over proteins...")
  res <- t(apply(what, 1, function(protein)
    {

      peptides <- as.numeric(unlist(protein["pep.ids"]))
      prot.id <- unlist(protein["prot.id"])

      all_peptide_fcs <- as.vector(random_fcs[peptides, ])
      all_peptide_fcs <- all_peptide_fcs[is.finite(all_peptide_fcs)]
      all_peptide_fcs <- sort(all_peptide_fcs)
      log2FC <- all_peptide_fcs[0.5*length(all_peptide_fcs) + 1]

      FC_SIZES <- floor(FC_LIMITS * length(all_peptide_fcs))
      fc_ranges <- sapply(FC_SIZES, function(size)
      {
        range_starts <- all_peptide_fcs[1:(length(all_peptide_fcs)-size + 1)]
        range_ends <- all_peptide_fcs[size:length(all_peptide_fcs)]
        ranges <- range_ends - range_starts

        best_idx <- which.min(ranges)

        r <- c(range_starts[best_idx], range_ends[best_idx])

        return(r)
      })

      range <- fc_ranges[, 1] # 50 percent interval

      loc_medians <- medians[peptides]

      smaller <- (loc_medians < 0 & loc_medians < range[1])
      taller <-  (loc_medians > 0 & loc_medians > range[2])
      to_correct <- smaller | taller

      if(sum(to_correct) > 0)
      {

        loc_fcs_shifted <- matrix(fcs[peptides, ], ncol=ncol(fcs), nrow=length(peptides))

        if(sum(smaller) > 0)
        {
          loc_fcs_shifted[smaller, ] <- loc_fcs_shifted[smaller, ] - range[1]
        }
        if(sum(taller) > 0)
        {
        loc_fcs_shifted[taller, ] <- loc_fcs_shifted[taller, ] - range[2]
        }

        loc_shifted_scores <- toZ(peptides[to_correct], matrix(loc_fcs_shifted[to_correct, ], nrow=sum(to_correct), ncol=ncol(fcs)), assignment, dists)


        loc_diffs <- abs(test_scores[peptides[to_correct], 1]) - abs(loc_shifted_scores)
        loc_s <- s[peptides[to_correct]]
        loc_bg_normed_score <- norm_background(loc_diffs,
                                                pep2validpair[peptides[to_correct]],
                                                loc_s,
                                                background_distributions,
                                                background_zeros)

        bg_normed_score[peptides[to_correct]] <- apply(cbind(abs(bg_normed_score[peptides[to_correct]]), abs(loc_bg_normed_score)), 1, min)

        bg_normed_score[peptides[to_correct]] <- bg_normed_score[peptides[to_correct]] * loc_s
      }

      prot.s <- sum(test_scores[peptides, 1])
      prot.sd <- sqrt(sum(pep2fullvar[peptides]))
      prot.p.val <- 2 * (1.0 - pnorm(abs(prot.s), sd=prot.sd))


      ss <- sum(bg_normed_score[peptides])
      p.val <- 2 * (1.0 - pnorm(abs(ss), sd=sqrt(length(peptides))))
      return(c(prot.id, ss, p.val, log2FC, prot.p.val, prot.sd, prot.s))
  }))
  colnames(res) <- c("prot.id", "score", "p.val", "log2FC", "prot.p.val", "prot.sd", "prot.s")
  res <- as.data.frame(res)
  res$p.adj <- p.adjust(res$p.val, method="BH")
  res$prot.p.adj <- p.adjust(res$prot.p.val, method="BH")
  rownames(res) <- res$prot.id
  for(i in 2:ncol(res))
  {
    res[, i] <- as.numeric(res[,i])
  }
  message("finished analysis")
  return(res)
}
