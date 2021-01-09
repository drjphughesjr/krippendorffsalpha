
######################################################################
#
# zzz.R
#
# Written by John Hughes <drjphughesjr@gmail.com>.
#
# Last Modified 05/27/20
# Licensed under the GNU General Public License version 2 (June, 1991)
#
######################################################################

#' @importFrom utils packageDescription
#' @importFrom stats quantile
#' @importFrom stats influence
#' @importFrom stats confint
#' @importFrom stats density
#' @importFrom graphics plot
#' @importFrom graphics abline
#' @importFrom graphics hist
#' @importFrom graphics lines

.onAttach = function(libname, pkgname)
{
    temp = packageDescription("krippendorffsalpha")
    msg = paste(temp$Package, ": ", temp$Title, "\n", "Version ", temp$Version,
                " created on ", temp$Date, ".\n", sep = "")
    msg = paste(msg, "copyright (c) 2020-2021, John Hughes\n",
                sep = "")
    msg = paste(msg, 'For citation information, type citation("krippendorffsalpha").\n', sep = "")
    msg = paste(msg, 'Type help(package = krippendorffsalpha) to get started.\n', sep = "")
    packageStartupMessage(msg)
}

