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



num_covar <-function(pair)
{
    N <- pair[1]
    M <- pair[2]
    return((N * M * (N + M - 2)) / 2)
}

num_var <-function(pair)
{
  N <- pair[1]
  M <- pair[2]
  return(N*M)
}

#' @export
generate_distribution <- function(pair, SD=NULL, N=100000, out.dir=NULL, with.Z=T)
{
    if(is.null(SD))
    {
      N_to_draw = N
      SD <- runif(N_to_draw) * 10 + 0.2
    }

    c1 <- matrix(rnorm(N * pair[1], 0, SD), nrow=N, ncol=pair[1])
    c2 <- matrix(rnorm(N * pair[2], 0, SD), nrow=N, ncol=pair[2])

    pairs <- t(expand.grid(1:pair[1], 1:pair[2]))

    SD_diff <- sqrt(2) * SD #since the fold changes are a linear combination of two normal factors

    #compute scores and converto to standar normal
    difference <- apply(pairs, 2, function(p) c1[, p[1]] - c2[, p[2]])
    if(with.Z==T)
    {
      #sum of Z values
      stats <- rowSums(difference)
      stats <- stats / SD_diff
    } else
    {
      stats <- pnorm(difference, sd=SD_diff) #transorfm to probability
      stats <- stats - 0.5 # shift such that values below median get a negative score, raning from 0 to -0.5
      stats <- rowSums(stats)
    }

    #apply factor to minimize differences
    factors <- log2(rowSums(2^c1)) - log2(rowSums(2 ^ c2)) - log2(pair[1]) + log2(pair[2])
    c2 <- c2 + factors

    #compute scores for shifted data
    sim_difference <- apply(pairs, 2, function(p) c1[, p[1]] - c2[, p[2]])
    
    if(with.Z==T)
    {
      #sum of Z values
      sim_stats <- rowSums(sim_difference)
      sim_stats <- sim_stats / SD_diff
    } else
    {
      sim_stats <- pnorm(sim_difference, sd=SD_diff)
      sim_stats <- sim_stats - 0.5
      sim_stats <- rowSums(sim_stats)
    }

    # this is the difference of two N(0,1) scores, one shifted, one not
    #compute differences and remove if shifted score is more extreme than real score
    stat_diffs <- abs(stats) - abs(sim_stats)
    stat_diffs[stat_diffs < 0 ] <- 0

    if(!is.null(out.dir))
    {
      pdf(file.path(out.dir, paste0(pair[1], ".vs.", pair[2], ".pdf")))
      plot(ecdf(stat_diffs))
      dev.off()
    }

    stat_diffs <- sort(stat_diffs)
    return(stat_diffs)
}
