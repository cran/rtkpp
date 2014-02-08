#-----------------------------------------------------------------------
#     Copyright (C) 2012-2014  Serge Iovleff, University Lille 1, Inria
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as
#    published by the Free Software Foundation; either version 2 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public
#    License along with this program; if not, write to the
#    Free Software Foundation, Inc.,
#    59 Temple Place,
#    Suite 330,
#    Boston, MA 02111-1307
#    USA
#
#    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
#
#-----------------------------------------------------------------------
#' rtkpp is a R/STK++ bridge
#'
#' This package contains the header files for the STK++ library.
#' The typical usage is to install this package and list it in
#' the \env{LinkingTo: } line in the \file{DESCRIPTION} file of
#' other packages.
#'
#' As described at the STK++ project's home page, \url{http://www.stkpp.org},
#' STK++ is a versatile, fast, reliable and elegant collection of C++ classes
#' for statistics, clustering, linear algebra, arrays (with an Eigen-like API),
#' regression, dimension reduction, etc. Some functionalities provided by the
#' library are available in the R environment as R functions.
#' 
#' The available functionalities are:
#' \enumerate{
#'    \item the clusterDiagGaussian method allowing to discover, group structures in
#' multivariate data sets using diagonal Gaussian multidimensional models (8 models).
#'    \item the clusterGamma method allowing to discover, group structures in
#' multivariate data sets using gamma multidimensional models (24 models).
#'  }
#'
#' \tabular{ll}{
#'   Package: \tab rtkpp\cr
#'   Type: \tab Package\cr
#'   Version: \tab 0.8.1\cr
#'   Date: \tab 2014-07-05\cr
#'   License: \tab GPL for the rtkpp side, LGPL for the stkpp side  + file LICENSE\cr
#'   LazyLoad: \tab yes\cr
#' }
#'
#' @rdname rtkpp-package
#' @name rtkpp
#' @aliases rtkpp
#' @docType package
#' @keywords STK++, stkpp
#' @import Rcpp
#'
#' @author
#' Author: Serge Iovleff \email{contact@@stkpp.org}
#'
#' @useDynLib rtkpp
NULL

#' Quantitative data: Old Faithful Geyser
#'
#' The file geyser.rda contains 272 observations from the Old Faithful Geyser
#' in the Yellowstone National Park. Each observation consists of two
#' measurements: the duration (in minutes) of the eruption and the waiting
#' time (in minutes) to the next eruption.
#'
#' Old Faithful erupts more frequently than any other big geyser, although it
#' is not the largest nor the most regular geyser in the park. Its average
#' interval between two eruptions is about 76 minutes, varying from
#' 45 - 110 minutes. An eruption lasts from 1.1/2 to 5 minutes,
#' expels 3,700 - 8,400 gallons (14,000 - 32,000 liters) of boiling water, and
#' reaches heights of 106 - 184 feet (30 - 55m). It was named for its consistent
#' performance by members of the Washburn Expedition in 1870. Old Faithful is
#' still as spectacular and predictable as it was a century ago.
#'
#' @format A data frame with 272 observations on the following 2 variables.
#'
#' \describe{
#'   \item{\code{Duration}}{a numeric vector containing the duration (in minutes) of the eruption}
#'   \item{\code{Waiting.Time}}{a numeric vector containing the waiting time (in minutes) to the next eruption}
#' }
#'
#' @source \url{http://www.geyserstudy.org/geyser.aspx?pGeyserNo=OLDFAITHFUL}
#'
#' @references
#' Hardle, W. (1991). "Smoothing Techniques with Implementation in S". Springer-Verlag, New York.
#' Azzalini, A. and Bowman, A. W. (1990). "A look at some data on the Old Faithful geyser". Applied Statistics 39, 357-365.
#'
#' @name geyser
#' @docType data
#' @keywords datasets
#'
#' @examples
#'   data(geyser)
NULL

#' Qualitative data: Survival of passengers on the Titanic
#' 
#' For each person on board the fatal maiden voyage of the ocean liner Titanic,
#' this dataset records: sex, age [adult/child], economic status [first/second/third class, or crew]
#' and whether or not that person survived. There are no missing values.
#'
#' The sinking of the Titanic is a famous event, and new books are still being
#' published about it. Many well-known facts-from the proportions of first-class
#' passengers to the "women and children first" policy, and the fact that that
#' policy was not entirely successful in saving the women and children in the
#' third class-are reflected in the survival rates for various classes of passenger.
#'
#' These data were originally collected by the British Board of Trade in their
#' investigation of the sinking. Note that there is not complete agreement among
#' primary sources as to the exact numbers on board, rescued, or lost.
#'
#' Due in particular to the very successful film "Titanic", the last years saw
#' a rise in public interest in the Titanic. Very detailed data about the
#' passengers is now available on the Internet, at sites such as
#' "Encyclopedia Titanica" (\url{http://www.rmplc.co.uk/eduweb/sites/phind}).
#'
#' @format A data frame with 2201 observations on the following 4 variables.
#'
#' \describe{
#'   \item{\code{Class}}{0 = crew, 1 = first, 2 = second, 3 = third, which denote the economic status of the subject}
#'   \item{\code{Age}}{1 = adult, 0 = child, which denote if the subject is an adult or a child}
#'   \item{\code{Sex}}{1 = male, 0 = female, which denote the sex of the subject}
#'   \item{\code{Survived}}{1 = yes, 0 = no, which denote if the subject lived through the fatal maiden voyage of the ocean liner Titanic} 
#' }
#'
#' @source 
#' The source provides a data set recording class, sex, age, and survival status
#' for each person on board of the Titanic, and is based on data originally
#' collected by the British Board of Trade and reprinted in:
#' British Board of Trade (1990), "Report on the Loss of the Titanic (S.S.)".
#'  British Board of Trade Inquiry Report (reprint).  Gloucester, UK: Allan Sutton Publishing.
#'
#' @name titanic
#' @docType data
#' @keywords datasets
#' 
#' @examples
#'   data(titanic)
NULL

#' Qualitative data : morphological description of birds
#' 
#' The dataset contains details on the morphology of birds (puffins). Each
#' individual (bird) is described by 6 qualitative variables. One variable for
#' the gender and 5 variables giving a morphological description of the birds.
#' There is 69 puffins divided in 2 sub-classes: lherminieri (34) and subalaris (35).
#'
#' @format A data frame with 69 observations on the following 5 variables.
#'
#' \describe{
#'   \item{\code{gender}}{a character vector defining the gender (2 modalities, male or female).}
#'   \item{\code{eyebrow}}{a character vector describing the eyebrow stripe (4 modalities).}
#'   \item{\code{collar}}{a character vector describing the collar (5 modalities).}
#'   \item{\code{sub-caudal}}{a character vector describing the sub-caudal (5 modalities).}
#'   \item{\code{border}}{a character vector describing the border (3 modalities).}
#' }
#'
#' @source
#' Bretagnolle, V., 2007. Personal communication, source: Museum.
#'
#' @name birds
#' @docType data
#' @keywords datasets
#' 
#' @examples
#'   data(birds)
NULL


#' Qualitative data : Car Evaluation
#' 
#' Car Evaluation Database was derived from a simple hierarchical decision model
#' originally developed for the demonstration of DEX, M. Bohanec, V. Rajkovic:
#' Expert system for decision making.
#'
#' @format A data frame with 1728 observations on the following 6 variables.
#'
#' \describe{
#'   \item{\code{buying}}{the buying price (4 modalities: vhigh, high, med, low).}
#'   \item{\code{maint}}{the price of the maintenance (4 modalities: vhigh, high, med, low).}
#'   \item{\code{doors}}{the number of doors (4 modalities: 2, 3, 4, 5more).}
#'   \item{\code{persons}}{the capacity in terms of persons to carry (3 modalities: 2, 4, more).}
#'   \item{\code{lug_boot}}{the size of luggage boot  (3 modalities: small, med, big).}
#'   \item{\code{safety}}{the estimated safety of the car (3 modalities: low, med, high).}
#'   \item{\code{acceptability}}{the car acceptability (4 modalities: unacc, acc, good, vgood).}
#' }
#'
#' @source
#' Creator: Marko Bohanec 
#' Donors: Marko Bohanec & Blaz Zupan
#' http://archive.ics.uci.edu/ml/datasets/Car+Evaluation
#'
#' @name car
#' @docType data
#' @keywords datasets
#' 
#' @examples
#'   data(car)
NULL

