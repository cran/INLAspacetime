% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/worldMap.R
\name{worldMap}
\alias{worldMap}
\title{Get the world map}
\usage{
worldMap(
  crs = sp::CRS("+proj=moll +units=km"),
  scale = "medium",
  returnclass = c("sp", "sf")
)
}
\arguments{
\item{crs}{a CRS object for the target coordinates.
Default is the Mollweide projection with units in kilometers.}

\item{scale}{see the help from \code{\link[rnaturalearth:ne_countries]{rnaturalearth::ne_countries()}}}

\item{returnclass}{defined by the \link{sp} if returnclass = "sp" or
by the \link{sf} package if returnclass = "sf".}
}
\value{
the world map
}
\description{
Consider the world map, fix the wrap around and convert to the
Mollweide projection in km distance (default, can be changed).
}
\details{
Requires the \code{rnaturalearth} package
}
\references{
the world map data is from the maps R package.
}
