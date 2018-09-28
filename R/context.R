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

get_mean_wo <- function(c)
{
  x <- mean(c[c!=0])
  return(x)
}

get_SD <- function(x)
{
  idx_to_take <- floor(0.341 * length(x)) + 1
  tmp <- sort(abs(x))
  return(tmp[idx_to_take])
}

get_intercond_foldchanges <- function(data)
{
  tmp <- as.vector(apply(combn(1:ncol(data), 2), 2, function(pair)
    {
      f <- log2(data[, pair[2]] / data[, pair[1]])
      return(f)
  }))
  tmp <- tmp[is.finite(tmp)]
  return(tmp)
}

#' @export
build_context <- function(data, c1, c2, stepsize=0.25, nsteps=NULL, out.dir=NULL)
{
  message("building context...")
  x <- exprs <-Biobase::exprs(data)

  #extract mean signal
  m1 <- apply(x[, c1], 1, FUN=function(r) get_mean_wo(r))
  m2 <- apply(x[, c2], 1, FUN=function(r) get_mean_wo(r))
  m <- apply(cbind(m1,m2), 1, FUN=min)

  row_order <- order(m)

  #reorder data by mean expression
  x <- x[row_order, ]
  m <- m[row_order]

  if(is.null(nsteps))
  {
    nsteps <- min(600, floor(nrow(exprs) / 20))
  }

  distributions <- list()
  counter <- 1
  pep2dist <- rep(-1, nrow(x))
  distances <- rep(-1, nrow(x))
  sdsss <- c()
  starts <- c()
  ends <- c()
  start_i <- c()
  end_i <- c()
  sizes <- c()
  stepsize <- floor(stepsize * nsteps)

  for(i in seq(nsteps, nrow(x) -stepsize, stepsize))
  {
    #update pep2dist assignments
    fi <- i-nsteps+1
    li <- min(i+nsteps, length(m))

    start <- m[fi]
    end <- m[li]
    starts <-c(starts, start)
    ends <- c(ends, end)
    start_i <- c(start_i, fi)
    end_i <- c(end_i, li)

    #reassign peptides to distribution
    sapply(fi:li, function(index)
      {
        s <- m[index]
        distance <- min(abs(start - s), abs(end - s))
        if(distances[index] <= distance)
        {
          distances[index] <<- distance
          pep2dist[row_order[index]] <<- counter
        }

    })

    #select subset of peptides and add inner cond +/- fold changes
    sub1 <- x[fi:li, c1]
    fold_changes1 <- get_intercond_foldchanges(sub1)
    fold_changes1 <- c(fold_changes1, -fold_changes1)

    sub2 <- x[fi:li, c2]
    fold_changes2 <- get_intercond_foldchanges(sub2)
    fold_changes2 <- c(fold_changes2, -fold_changes2)


    # add fold change to mean to the changes
    sub_m1 <- apply(sub1, 1, FUN=function(r) get_mean_wo(r))
    sub_m2 <- apply(sub2, 1, FUN=function(r) get_mean_wo(r))


    tmp <- log2(sub1 / sub_m1)
    tmp <- tmp[is.finite(tmp)]
    fold_changes1 <- c(fold_changes1, tmp)

    tmp <- log2(sub2 / sub_m2)
    tmp <- tmp[is.finite(tmp)]
    fold_changes2 <- c(fold_changes2, tmp)

    # normalize and compare between conditions
    reverse <- sub_m2 < sub_m1
    sub1[reverse] <- sub1[reverse] * (sub_m2[reverse] / sub_m1[reverse])
    sub2[!reverse] <- sub2[!reverse] * (sub_m1[!reverse] / sub_m2[!reverse])


    fold_changes3 <- as.vector(apply(t(expand.grid(1:sum(c1), 1:sum(c2))), 2, function(pair)
      {
        tmp <- log2(sub2[, pair[2]] / sub1[, pair[1]])
        return(tmp)
      }))
    fold_changes3 <- fold_changes3[is.finite(fold_changes3)]

    final_distrib <- c()

    sd1n <- get_SD(fold_changes1[fold_changes1 <= 0])
    sd2n <- get_SD(fold_changes2[fold_changes2 <= 0])
    sd3n <- get_SD(fold_changes3[fold_changes3 <= 0])
    max_sd <- order(c(sd1n, sd2n, sd3n))[3]
    if(max_sd==1)
    {
      final_distrib <- c(final_distrib, fold_changes1[fold_changes1 <= 0])
    }
    else if(max_sd==2)
    {
      final_distrib <- c(final_distrib, fold_changes2[fold_changes2 <=0])
    }
    else
    {
      final_distrib <- c(final_distrib, fold_changes3[fold_changes3 <=0])
    }

    sd1p <- get_SD(fold_changes1[fold_changes1 >=0])
    sd2p <- get_SD(fold_changes2[fold_changes2 >=0])
    sd3p <- get_SD(fold_changes3[fold_changes3 >=0])
    max_sd <- order(c(sd1p, sd2p, sd3p))[3]
    if(max_sd==1)
    {
      final_distrib <- c(final_distrib, fold_changes1[fold_changes1 >=0])
    }
    else if(max_sd==2)
    {
      final_distrib <- c(final_distrib, fold_changes2[fold_changes2 >=0])
    }
    else
    {
      final_distrib <- c(final_distrib, fold_changes3[fold_changes3 >=0])
    }
    final_distrib <- sort(final_distrib)

    if(!is.null(out.dir))
    {
      pdf(file.path(out.dir, paste0(counter, "_distrib.pdf")))
      plot(density(final_distrib), main=sprintf("size: %d", length(final_distrib)))
      dev.off()
    }

    sizes <- c(sizes, length(final_distrib))
    distributions[[counter]] <- final_distrib

    idx_to_take <- floor(length(final_distrib) * 0.682) + 1
    sdsss <- c(sdsss, sort(abs(final_distrib))[idx_to_take])

    counter <- counter + 1
  }
  pep2dist[pep2dist == -1] <- length(distributions)
  message(sprintf("created %d context distributions", length(distributions)))
  return(list(distributions=distributions, assignments=pep2dist, sds=sdsss, starts=starts, ends=ends,
              start_i=start_i, end_i=end_i, sizes=sizes))
}
