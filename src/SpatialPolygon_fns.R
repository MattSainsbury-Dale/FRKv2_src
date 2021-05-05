## Lower left-hand vertex is (a, b), side length s
makeSquare <- function(a, b, s) {
  
  ## Start with lower left-hand vertex then continue clockwise
  vertex_corners <- rbind(c(a, b), 
                          c(a, b + s), 
                          c(a + s, b + s), 
                          c(a + s, b))
  
  square <- Polygon(vertex_corners)
  
  return(square)
}

## s is the side length of the grid squares, 
## n is the number of grid squares in each dimension,
makeGrid <- function (n, range0 = 1) { 
  
  ## range0 is the desired extent of the grid.
  ## Define s, the side length of the squares, so that s*n = range = range0:
  s <- range0 / n
  
  squares <- list()
  i <- 1
  for (b in seq(0, (n-1)*s, by = s)) {
    for (a in seq(0, (n-1)*s, by = s)) {
      square <- Polygons(list(makeSquare(a = a, b = b, s = s)), 
                         paste0("square", i))
      squares[[i]] <- square
      i = i + 1
    }
  }
  
  square_grid <- SpatialPolygons(squares)
  
  return(square_grid)
}


## Shift coordinates by a: a is scalar or length-two vector
shift.coordinates <- function(coords, a) {
  
  if (length(a) == 1) {
    coords <- coords + a
  } else if (length(a) == 2) {
    coords <- coords + rep(a, each = nrow(coords))
  } else 
    stop("a should be a scalar or a vector of length 2.")
}


## Credit to this function goes to:
# https://stackoverflow.com/a/4860337
convex.poly <- function(nSides, area) 
{
  # Find the radius of the circumscribed circle, and the angle of each point if this was a regular polygon
  radius <- sqrt((2*area)/(nSides*sin((2*pi)/nSides)))
  angle <- (2*pi)/nSides
  
  # Randomize the radii/angles
  radii <- rnorm(nSides, radius, radius/10)
  angles <- rnorm(nSides, angle, angle/10) * 1:nSides
  angles <- sort(angles)
  
  points <- list(x=NULL, y=NULL)
  points$x <- cos(angles) * radii
  points$y <- sin(angles) * radii
  
  # Find the area of the polygon
  m <- matrix(unlist(points), ncol=2)
  m <- rbind(m, m[1,])
  current.area <- 0.5 * (sum(m[1:nSides,1]*m[2:(nSides+1),2]) - sum(m[1:nSides,2]*m[2:(nSides+1),1]))
  
  points$x <- points$x * sqrt(area/current.area)
  points$y <- points$y * sqrt(area/current.area)
  
  ## Append the first set of points to the end, to "complete" the polygon
  points$x <- append(points$x, points$x[1])
  points$y <- append(points$y, points$y[1])
  
  ## Convert to matrix
  points <- do.call(cbind, points)
  
  return (points)
}

## coords: list of matrices containing the coordinates of the polygons
## This function is based on the simple example from vignette("sp")
coords_to_SpatialPolygons <- function(coords) {
  
  ## Convert to Polygon
  Poly_list <- lapply(coords, Polygon)
  
  ## Give each Polygon a name
  names(Poly_list) <- paste0("s", 1:length(Poly_list))
  
  ## Convert to Polygons
  Polygons_list <- lapply(seq_along(Poly_list), function(i) Polygons(list(Poly_list[[i]]), names(Poly_list)[i]))
  
  ## Convert to SpatialPolygons
  newdata <- SpatialPolygons(Polygons_list, 1:length(Polygons_list))
  
  return(newdata)
}