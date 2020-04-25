
"
GeomScale Medium Challenge GSOC 2020

Implementation of the Randomised Cuting Plane(RCP) algorithm 

Referenced paper:
Dabbene, Fabrizio, Pavel S. Shcherbakov, and Boris T. Polyak. 
A randomized cutting plane method with probabilistic geometric convergence 
SIAM Journal on Optimization 20.6 (2010): 3185-3207.
"


library(volesti)
library(Rglpk)
library(lpSolve)
library(truncnorm)
library(ggplot2)



RCP_algorithm <- function(A,b,c,epsilon=1e-8,n_loops=100,accuracy=1e-8){
  "
  Function for solving linear programming problems, i.e:
                        min{t(c)*x}
                        Ax<=b
  where t(c) is the transpose matrix of c.
  The condition x>=0 must be inserted in matrices A,b.

  Parameters:
  A: matrix with left constraints for x
  b: vector of dimension nx1 with right constraints for x
  c: vector of dimesion 1xn included in function to be minimized
  epsilon: parameter for number of samples at each iteration
  n_loops: max number of iterations
  accuracy: accuracy/tolerance of result

  Output:the optimum of the function t(c)*x
  
  "
  
  #Initializing variables
  min_val <- 0 #initial optimum
  polytope <- Hpolytope$new(A,b)  
  k <- 1  #number of iteration
  error <- 100 #difference between previous optimum and current optimum
  x_cur <- matrix(rtruncnorm(n=length(A[1,]),a=0),length(A[1,]),1) #random initial point x_k
  flag <- 0 #activates if x<0 is generated signaling to stop
  
  
  while(k < n_loops && error > accuracy){
    
    #find num_samples for this iteration
    num_samples <- ceiling( 2.2 * log(1/epsilon) + 1.1 + 0.505 * log(k) )
    
    
    #sample from polytope
    points <- sample_points(P=polytope,N=num_samples,WalkType = "CDHR")
    
    
    #save previous values
    x_old <- x_cur
    prev_min <- min_val
    
    #find optimum value from generated samples
    index <- which.min(t(c) %*% points)
    x_cur <- points[,index]
    min_val <- t(c) %*% x_cur
    
    #restrict polytope 
    polytope$A <- rbind(polytope$A, t(c))
    polytope$b <- c(polytope$b, min_val)
    
    
    #Find error between prev_min and cur_min
    error <- abs(as.vector(min_val) - as.vector(prev_min))
    
    k <- k + 1
    
    if(all(is.infinite(points) || is.na(points))) {
      min_val <- prev_min
      x_cur <- x_old
      break
    }
    
  }
  #if flag activated make x=0
  for(i in 1:length(x_cur)){
    if(x_cur[i]<0)
      flag <-1
  }
  if(flag) x_cur=matrix(0,length(x_cur),1)
  
  
  return(min_val)
}



test_fun <- function(){
  "
  Function for running tests making benchmarks of the RCP algorithm implementation.

  Outputs: a csv file with benchmarks
  "
  
  #Initialising parameters
  result_error <- "Error"
  rcp_ans <- "RCP Algorithm"
  lp_ans <- "lpSolve() Algorithm"
  times_lp <- "lpSolve Time"
  times_rcp <- "RCP time"
  
  #number of digits to be returned by sys.time()
  op <- options(digits.secs = 6) 
  options(op)
  
  dimensions = c(rep(5,10),rep(10,10),rep(20,10))
  
  for(i in 1:30){
    
    current_dim <- dimensions[i]
    
    #Create random matrices A,b,c
    A <- matrix(runif(n=current_dim*current_dim,min=-20,max=20),current_dim,current_dim)
    b <- matrix(runif(n=current_dim,min=0,max=20),current_dim,1)
    c <- matrix(runif(n=current_dim,min=-20,max=20),current_dim,1)
    Amat <- A
    bmat <- b
    
    #Add constraint x>=0
    A <- rbind(A, diag(-1,current_dim,current_dim))
    b <- c(b, matrix(0,current_dim,1))
    
    #Find solution with lpSolve() and time needed for it
    start_time <- Sys.time()
    
    lp_optimum <- checker(Amat,bmat,c)
    lp_optimum <- signif(lp_optimum,digits=6)
    
    end_time <- Sys.time()
    
    #Find time needed for lpSolve() result
    time <- end_time - start_time
    time <- signif(time,digits=6)
    times_lp <- c(times_lp, time)
    
    
    lp_ans <- rbind(lp_ans,lp_optimum)
    
    #Find solution with RCP algorithm and time needed for it
    start_time <- Sys.time()
    rcp_optimum <- RCP_algorithm(A,b,c)
    rcp_optimum<- signif(rcp_optimum,digits=6)
    end_time <- Sys.time()
    
    #Find time needed for RCP algorithm result
    time <- end_time - start_time
    time <- signif(time,digits=6)
    times_rcp <- c(times_rcp, time)
    
    
    rcp_ans <- rbind(rcp_ans,rcp_optimum)
    
    #Find difference between results from RCP and lpSolve()
    res_error <- abs(abs(lp_optimum) - abs(rcp_optimum))
    res_error <- signif(res_error,digits=6)
    result_error <- rbind(result_error,res_error)
  }
  
  
  time_arr <- cbind(times_lp,times_rcp)
  result <- cbind(lp_ans,rcp_ans,result_error,times_lp,times_rcp)
  
  
  #create csv file with benchmarks
  write.csv(result, file = "benchmarks_medium_task2.csv", row.names = FALSE)
  
  
  #plot time needed for RCP vs lpSolve() 
  g <- ggplot(data.frame( x=1:31, y=times_lp[1:length(times_lp)] ),aes(x="lpSolve()", y=y, colour="lpSolve() time")) +
    geom_point( aes(x=x, y=y, color = "lpSolve() time")) +
    coord_fixed(xlim = c(1,31)) +
    ggtitle("Time of RCP vs lpSolve()")+geom_point(aes(x,y = times_rcp[1:length(times_lp)], color="RCP time"))+ xlab("dimensions")+ylab("Time")
    
  plot(g)
  
}

checker <- function(A,b,c){
  "
  Function for solving a min linear program with lpSolve().
  
  Parameters:
  A: left array constraints(without the constraints for x>=0)
  b: right array constraints(without the constraints for x>=0)
  c: array included in function to be minimized
  
  Outputs: the optimum value of function
  
  Reference: https://towardsdatascience.com/linear-programming-in-r-444e9c199280
  "
  
  # Set coefficients of the objective function
  f.obj <- c
  
  # Set matrix corresponding to coefficients of constraints by rows
  # Do not consider the non-negative constraint; it is automatically assumed
  f.con <- A
  
  # Set unequality signs
  f.dir <- c(replicate(dim(b)[1],"<="))
  
  # Set right hand side coefficients
  f.rhs <- b
  
  # Final value (z)
  #optimum <- lp("min", f.obj, f.con, f.dir, f.rhs)
  z <- lp("min", f.obj, f.con, f.dir, f.rhs)$solution
  optimum <- t(c) %*% z
  
  return(optimum)
}

test_fun()
