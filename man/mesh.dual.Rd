% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mesh.dual.R
\name{mesh.dual}
\alias{mesh.dual}
\title{Extracts the dual of a mesh object.}
\usage{
mesh.dual(mesh, SP = TRUE, mc.cores = getOption("mc.cores", 2L))
}
\arguments{
\item{mesh}{a 2d mesh object.}

\item{SP}{logical indicating if the output is to be returned
as a SpatialPolygons object. Default is TRUE.}

\item{mc.cores}{number of threads to be used.}
}
\value{
a list of polygons or a SpatialPolygons object.
}
\description{
Extracts the dual of a mesh object.
}
