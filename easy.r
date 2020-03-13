# install.packages('volesti')
# install.packages('ggplot2')

library(ggplot2)
library(volesti)

# sample the 2-dimensional cube

P <- GenCube(2, 'H')

for (walk in c("CDHR", "RDHR", "BW")) {
  
  points <- sample_points(P, WalkType = walk, N = 100)
  
  g <- ggplot(data.frame( x=points1[1,], y=points1[2,] )) + 
              geom_point(aes(x=x, y=y, color=walk))
  plot(g)
}
