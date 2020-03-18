  library(volesti)
  library(ggplot2)
  
  P <- GenCube(100, 'H')
  points <- sample_points(P,N=1000)
  g <- ggplot(data.frame( x=points[1,], y=points[2,] )) +
              geom_point( aes(x=x, y=y, color="red")) + 
              coord_fixed(xlim = c(-1,1),ylim = c(-1,1)) + 
              ggtitle("Sampling from 10 dimensional hypercube")
  plot(g)
  