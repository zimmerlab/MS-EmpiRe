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


get_error <- function(mat, c1=NULL, c2=NULL)
{
    fcs <- generate_foldchanges(mat, c1, c2, remove_na=TRUE, as_v=TRUE)
    return(sd(fcs))
}


shift <- function(mat, c1 , c2)
{
    shift_par <- median(generate_foldchanges(mat, c1, c2, remove_na=TRUE, as_v=TRUE))
    return(mat[, c2] * (2 ^ (-shift_par)))
}

cluster <- function(mat)
{
    message("starting clustering...")
    clusters <- 1:ncol(mat)
    for(i in 1:(ncol(mat) - 1))
    {
        ind_clusters <- unique(clusters)
        pairs <- combn(ind_clusters, 2)

        #compute distances
        distances <- apply(pairs, 2, function(pair) get_error(mat, clusters==pair[1], clusters==pair[2]))

        #select two nearest
        ordered <- order(abs(distances))
        p1 <- pairs[,ordered[1]]
        #p1 <- pairs[,which.min(abs(distances))]

        cluster1 <- p1[1]
        cluster2 <- p1[2]

        cluster2_size <- sum(clusters==cluster2)
        mat[, clusters==cluster2] <- shift(mat, clusters==cluster1, clusters==cluster2)

        clusters[clusters==cluster2] <- rep(cluster1, cluster2_size)
    }
    print(sprintf("final error: %g", get_error(mat)))
    return(mat)
}

plot_foldchanges <- function(data, out.file, c1=NULL, c2=NULL, mode=NULL)
{
    pdf(out.file)
    plot(c(),c(), type='n', xlim=c(-3,3), ylim=c(0, 1))

    pairs <- combn(1:ncol(data), 2)

    t1 <- data
    t2 <- data

    if(!is.null(c1))
    {
        t1 <- data[, c1]
        t2 <- data[, c2]

        pairs <- t(expand.grid(1:sum(c1), 1:sum(c2)))
    }


    colors <- 1:length(pairs)

    i <- 1
    apply(pairs, 2, function(pair)
      {
      fcs <- log2(t2[, pair[2]] / t1[, pair[1]])
      fcs <- fcs[is.finite(fcs)]
      lines(ecdf(fcs), col=colors[i])
      i <<- i+ 1
    })

    legend("bottomright", legend=apply(pairs, 2, function(pair) sprintf("%d vs %d", pair[1], pair[2])), col=colors)
    abline(v=0, lty=2)
    if(!is.null(mode))
    {
      abline(v=mode, lty=2)
    }
    abline(h=0.5, lty=2)

    dev.off()
}

detect_mode <- function(inputs, fc_width)
{
    x <- sort(inputs)
    half_width <- 0.5 * fc_width

    w_s <- t(apply(cbind(x[2:(length(x)-1)], 2:(length(x)-1)), 1, function(r)
      {
      fc<- r[1]
      idx <- r[2]

      start <- max(1, bin_search(x, fc-half_width, 1, idx-1) + 1)
      end <- min(length(x), bin_search(x, fc+half_width, idx + 1, length(x)) + 1)

      return(c(end-start, fc))
    }))

    return(w_s[which.max(w_s[,1]), 2])
}

#' @export
normalize <- function(data, out.dir=NULL)
{
  message("starting normalization...")
  conditions <- extract_conditions(data)
  cond1 <- conditions[,1] == 1
  cond2 <- conditions[,2] == 1

  x <- exprs(data)

  if(!is.null(out.dir))
  {
    message("plotting foldchanges pre normalization...")
    plot_foldchanges(x[, cond1], file.path(out.dir, "cond1.unnormed.pdf"))
    plot_foldchanges(x[, cond2], file.path(out.dir, "cond2.unnormed.pdf"))
    plot_foldchanges(x, file.path(out.dir, "combined.unnormed.pdf"), cond1, cond2)
  }

  message("starting to cluster...")
  x[, cond1] <- cluster(x[, cond1])
  x[, cond2] <- cluster(x[, cond2])



  SD <- max(get_error(x[, cond1]), max(get_error(x[, cond2])))
  print(get_error(x[, cond1]))
  print(get_error(x[, cond2]))
  fcs <- generate_foldchanges(x, cond1, cond2, remove_na=TRUE, as_v=TRUE)
  print("detecting mode")
  mode <- detect_mode(fcs, SD)
  print(sprintf("mode: %.g", mode))

  if(!is.null(out.dir))
  {
    message("plotting foldchanges after condition normalization")
    plot_foldchanges(x[, cond1], file.path(out.dir, "cond1.intranormed.pdf"))
    plot_foldchanges(x[, cond2], file.path(out.dir, "cond2.intranormed.pdf"))
    plot_foldchanges(x, file.path(out.dir, "combined.intranormed.pdf"), cond1, cond2, mode=mode)
  }



  x[, cond2] <- x[, cond2] * (2 ^ (-mode))

  if(!is.null(out.dir))
  {
    print(cond1)
    print(cond2)
    plot_foldchanges(x[, cond1], file.path(out.dir, "cond1.normed.pdf"))
    plot_foldchanges(x[, cond2], file.path(out.dir, "cond2.normed.pdf"))
    plot_foldchanges(x, file.path(out.dir, "combined.normed.pdf"), cond1, cond2, mode=mode)
  }
  print(sprintf("Final error: %g", get_error(x, cond1, cond2)))
  exprs(data) <- x
  return(data)
}
