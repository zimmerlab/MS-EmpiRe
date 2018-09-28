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



#' Binary search helper function
#' use bin_search instead
.bin_search_rec <- function(l, val, start, end)
{
    s <- start
    e <- end

    while(s <= e)
    {
        m <- bitwShiftR(s + e, 1)
        if(l[m] < val)
        {
          s <- m + 1
          next
        }
        if(l[m] > val)
        {
            e <- m - 1
            next
        }
        return(m + 1)
    }
    return(s)
}

#' Find a value in a list using binary search
#' Wrapper for bin_search_rec
#'
#' @param l A list of values
#' @param val A value that is looked for
#' @return Postition - 1 at which \code{val} must be in inserted to keep \code{l} sorted. \code{0} if \code{val} is
#' smaller than the first element of \code{l}, \code{length(l)} if it is greater than the last one
#' @examples
#' bin_search(1:10, 2)
#' bin_search(1:10, 0)
#' bin_search(1:10, 2.5)
bin_search <- function(l, val, start=NULL, end=NULL)
{
    if(is.infinite(val) || is.na(val) || is.null(val))
    {
      return(NA)
    }
    # if(val < l[1])
    # {
    #   return(1)
    # }
    #
    # if(val > l[length(l)])
    # {
    #   return(length(l) + 1)
    # }
    if(is.null(start))
    {
      start <- 1
    }
  if(is.null(end))
  {
    end <- length(l)
  }

    return(.bin_search_rec(l, val, start, end))
}

# convert pData to DesignMatrix
#' @export
extract_conditions <- function(data)
{
  condition <- factor(pData(data)$condition)
  return(model.matrix(~condition + 0))
}


nw_factor <- 4.0
n <- floor(1/0.05)
LIMITS <- sapply(1:(n-1), function(d)
{
  1.0 - (d / n)
})
weights <- dnorm(qnorm(LIMITS, sd=2), sd=2)
m <- min(weights)
iweights <- sapply(weights, function(d) round(nw_factor * d / m))
sampling <- matrix(c(LIMITS, iweights), ncol=2, nrow=length(LIMITS))

to_normal <- function(fc, sd)
{
  r <- as.numeric(as.vector(sapply(fc, function(f)
    {
      if(!is.finite(f))
      {
        return(rep(NaN, sum(sampling[,2])))
      }
      unlist(apply(sampling, 1, function(ff)
        {
          times <- ff[2]
          l <- ff[1]
          return(rep(qnorm(l, mean=f, sd=sd), times))
        }))

    })))
  return(r)
}
