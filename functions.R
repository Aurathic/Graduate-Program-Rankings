#!usr/bin/R

### This contains the basic functions that can be used to rank schools.

### These ones make and clean (optionally) a cross table of appointments

make.appointments <- function(data,min.year = 1990,max.year = 2005) {
  #This function requires that you have variables in the dataframe
  # named 'Year', 'School', and 'DegreeSchool' for the appropriate
  #variables
  sample <- data[(data$Year >= min.year),]
  sample <- sample[(sample$Year <= max.year),]
  appointments <- table(sample$DegreeSchool, sample$School)
  ###  Trim appointments list to schools where joint information is available
  appointments <- appointments[,colSums(appointments) >= 1]
  appointments <- make.appointments.symmetrical(appointments)
  appointments
}

make.appointments.symmetrical <- function(Appointments) {
  names	   	 <- intersect(dimnames(Appointments)[[1]], dimnames(Appointments)[[2]])
  Appointments <- Appointments[dimnames(Appointments)[[1]] %in% as.matrix(names),
                               dimnames(Appointments)[[2]] %in% as.matrix(names)]
  Appointments
}


###  regular.stochastic function makes a matrix regular stochastic by eliminating
###  any rows which sum to 0. It automatically
###  will stabilize within the same number of iterations as there are rows.

regular.stochastic  <- function (matrix) {
  names  <- dimnames(matrix)[[1]]
  for (i in 1:nrow(matrix)) {
    matrix <- matrix[dimnames(matrix)[[1]] %in% as.matrix(names),
                     dimnames(matrix)[[2]] %in% as.matrix(names)]
    names	 <- names[rowSums(matrix) > 0]
    names      <- as.matrix(names)
  }
  matrix
}

###  scale.rows function scales all rows against the number of degrees awarded
###  by a program.
scale.rows   <- function (Appointments,degrees) {
  degrees <- degrees[names(degrees) %in% dimnames(Appointments)[[1]]]
  for (i in 1:nrow(Appointments))  {
    
    Appointments[i,] <- Appointments[i,]/degrees[i]
  }
  Appointments
}

###  norm.columns function scales all columns so that each program begins with
###  one vote.
norm.columns <- function(Appointments) {
  for (i in 1:ncol(Appointments))  {
    if (sum(Appointments[,i]) >=1) {
      Appointments[,i] <- Appointments[,i]/sum(Appointments[,i],na.rm=T)
    }
  }
  Appointments
}
###  pagerank function determines relative influence of schools
pagerank <- function(Appointments, q=.1, weighting = "schools",
                     regular.stochastic=F, row.scaled=F,row.scalar=degrees) {
  
  if (regular.stochastic==T) {Appointments <- regular.stochastic(Appointments)}
  if (row.scaled==T)  {Appointments <- scale.rows(Appointments,degrees)}
  Appointments <- norm.columns(Appointments)
  #	Appointments <- norm.columns(Appointments)
  
  N   <- ncol(Appointments)     ## number of schools
  R   <- matrix(.1,N,1)	    ## Scores of the schools. Starts equal.
  rownames(R)  <- dimnames(Appointments)[[1]] ### Gives names to R
  
  
  ### Matrix one is the transitional matrix on a random selection.
  ### Weighting determines how the probabilites are determined.
  
  if (weighting=="degrees") {
    degrees2  <- degrees[names(degrees) %in% dimnames(Appointments)[[2]]]
    matrix1 <- matrix(degrees2/sum(degrees2)*q,N,1)}
  if (weighting=="faculty") {
    positions2<- positions[names(positions) %in% dimnames(Appointments)[[2]]]
    matrix1 <- matrix(positions2/sum(positions2)*q,N,1)}
  if (weighting=="schools") {matrix1 <- matrix(rep(q/N,N),N,1)}
  
  ### This is the actual formula that does the work. Applied recursively
  ### until the matrices stabilize.
  oldR         <- R+5  ### Keeps the matrix from stabilizing immediately
  
  while (abs(sum(oldR-R)) > .0000000000001) {
    oldR <- R
    R    <- (matrix1 + (((1-q) * Appointments) %*% R))
  }
  # Re-imagine as probabilities.
  R = R/sum(R)
  # n = n + 1; oldR <- R; R<-(matrix1 + (((1-q) * Appointments) %*% R));cbind(rownames(as.data.frame(R))[orderR)],sort(R))
  R
}
###  eigenrank returns the scores based on solving for the steady-state eigenvector
eigenrank <- function(matrix,row.scaled=T) {
  matrix  		  <- regular.stochastic(matrix)
  if (row.scaled==T)  {matrix <- scale.rows(matrix,degrees)}
  matrix  		  <- norm.columns(matrix)
  names   		  <- dimnames(matrix)[[1]]
  eigenvalues      <- eigen(matrix)[[1]]
  R 			  <- Re(eigen(matrix)$vectors[,which(round(Re(eigenvalues),15)==1)])
  names(R)         <- names
  R                <- as.matrix(R)
  dimnames(R)[[1]] <- names
  R
}

### Quick function to transform list of numbers to a scaled list with the top
### score normed to 100

pointscale <- function(matrix,degreescale=degrees, pointscale.cutoff = cutoff) {
  names    		<- rownames(matrix)
  localdegrees   <- degreescale[names(degreescale) %in% names]
  vector   		<- as.numeric(matrix[,1]) *100
  max 			<- max(vector[localdegrees >= pointscale.cutoff])
  min      		<- min(vector[localdegrees>= pointscale.cutoff])
  score    		<- 100 * (vector-min)/(max-min)
  names(score) 	<- names
  score
}

