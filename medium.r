Calc_number_of_samples <- function(k, eps) {
  #
  #
  # Implements the computation of the number of samples required at step k, with parameter eps following [B1, Section 4.2.2]
  #
  # Parameters:
  #
  # k: integer representing number of iteration.
  # eps: probability error tolerated.
  #
  # References:
  # [B1] Dabbene, Fabrizio, Pavel S. Shcherbakov and Boris T. Polyak. "A randomized cutting plane method with probabilistic geometric
  #      convergence". SIAM Journal of Optimization 20.6 (2010): 3185-3207.
  #

  ans <- 2.2 * log(1/eps) + 1.1 + 0.505 * log(k)
  
  return (ans)

}

Randomized_cutting_plane <- function(P, c, iter = 10000, eps = 1e-8, tol = 1e-8, positive = TRUE) {
  #
  #
  # Implements the randomized cutting algorithm from [B1] to solve
  # a linear programming problem of the form
  #
  #            min c * x 
  #            s.t. Ax <= b.   (1)
  #
  #
  # Parameters:
  # P: an H-polytope. The matrices A and b from the LP problem are, respectively, P$A and P$b.
  # c: row vector defining the objective function.
  # iter: the maximum number of iterations allowed to optimize. Its default value is iter = 1000.
  # eps: the probability involved in the probabilistic convergence rate of the algorithm (see [B1,Section 4.2]). Default value is eps = 1e-8.
  # tol: the tolerance for the minimum allowed improvement of the iterative solutions. Default value is tol = 1e-8
  # positive: a boolean to indicate if the variables are non-negative. The default value is positive = TRUE.
  #
  # Returns:
  #
  # ans: an object with two fields, solution and optimum. Solution is the vector which minimizes (1), while optimum is the value of the
  #      objective function at the vector solution.
  #
  # References:
  # [B1] [B1] Dabbene, Fabrizio, Pavel S. Shcherbakov and Boris T. Polyak. "A randomized cutting plane method with probabilistic geometric
  #      convergence". SIAM Journal of Optimization 20.6 (2010): 3185-3207.
  #

  # initialization
  
  k <- 1 # iteration
  sol <- NULL # current solution
  best <- NULL # current optimum
  stop <- FALSE # stopping condition
  
  # is positive is set to TRUE, to the initial polytope P0 we add the usual restrictions x >= 0
  
  dim <- P$dimension
  
  if(positive) {
    D <- diag(-1,dim,dim)
    P$A <- rbind(P$A, D)
    P$b <- c(P$b, matrix(0,dim,1))
  }
  
  # main loop  
  
  while(k <= iter && !stop) {
    
    # generate random samples from the current polytope

    N_k <- Calc_number_of_samples(k, eps) # calculates the number of samples N_k as in [B1, page]
    candidates <- sample_points(P, WalkType = "CDHR", N = N_k) 

	  # evaluate best solution among the sample points

    values <- c %*% candidates
    min_index <- which.min(values)
    arg_min <- candidates[,min_index]
    min_value <- c %*% arg_min

	  # update the polytope adding the constraint of the new cutting plane	

    P$A <- rbind(P$A, c)
    P$b <- c(P$b, min_value)
    
    # update current solution and optimum value
 
    sol <- arg_min
    prev_best <- best
    best <- min_value[1,1]

    # check stopping condition
    
    if( all(is.na(candidates) || is.infinite(candidates)) ){
      print("caca")
    }
    
    if(all(is.na(candidates) || is.infinite(candidates)) || 
           (!is.null(prev_best) && (prev_best - best) < tol) ) {
      stop <- TRUE
    }
    
    k <- k + 1
  }
  
  ans = list(solution = sol, optimum = best)
  
  return (ans)

}

Generate_tests<-function(filename) {
  #
  # 
  # Generates random linear programming problems using volesti::GenRandHpoly and
  # prints the results and running times to a file. In the file it also stores
  # the results obtained via the LP-solver in library Rglpk
  #
  # Parameters:
  # filename: name of the output file (should be a .csv file)
  #
  
  # initializations
  
  tests <- 5
  
  dimensions <- seq(10,50,10)
  facets <- seq(50,170,30)
  
  times_glpk <- c() # vector of executing times for GNU Linear Programming Kit, in seconds
  opt_glpk <- c() # vector of optima found by GNU Linear Programming Kit solver
  times_cut <- c() # vector of executing times for random cutting planes, in seconds
  opt_cut <- c() # vector of optima found by random cutting planes algorithm
  
  for (i in 1:5) {
    
    dim <- dimensions[i]
    fac <- facets[i]
    
    for (j in 1:tests) {
      
      # generate random LP problems. 
      
      P0 <- GenRandHpoly(dim, fac)
      obj <- runif(dim, -30, 30)
      
      # run GNU Linear Programming Kit (GLPK) solver on the problem
      
      dir <- rep("<=", fac) # direction of constraints; all of the type "<=".
      
      time <- system.time(x_glpk <- Rglpk_solve_LP(obj, P0$A, dir, P0$b))
      
      # update times and optima for GNU Linear Programming Kit solver
      
      times_glpk <- c(times_glpk, time)
      opt_glpk <- c(opt_glpk, x_glpk$optimum)
      
      # run randomized cutting algorithm implementation
      
      time <- system.time(x_cut <- Randomized_cutting_plane(P0, obj))
      
      # update times and optima for random cutting algorithm
      
      times_cut <- c(times_cut, time)
      opt_cut <- c(opt_cut, x_cut$optimum)
    }
  }
  
  # form data frame with all the data
  
  df <- data.frame(Dimension = rep(dimensions, each = tests), Facets = rep(facets, each = 5),
                     Time_GLPK = times_glpk, Optimum_GLPK = opt_glpk, 
                     Time_Cutting_Planes = times_cut, Optimum_Cutting_Planes = opt_cut)
  
  # print data to the file
  
  write.csv(df, filename, row.names = FALSE)
  
}

# install.packages('volesti')
# install.packages('Rglpk')

library('volesti')
library('Rglpk')

Generate_tests("~/data.csv")