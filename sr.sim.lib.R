#make.spines generates a spherical volume with given radius populated with synapses with given synapse density
make.spines <- function(synapse.density=2.07, radius=10){
  spine.radius <- 0.224/2 #http://dx.doi.org/10.1016/S0028-3908(98)00023-9 radius of average post synaptic density
  volume <- 4/3*pi*radius^3 #simulation volume
  n.synapses <- synapse.density * volume
  if(n.synapses<2){
    stop("Volume contains less than two spines. Try increasing the simulation volume or synapse.density.")
  }
  synapse.positions <- array(0,dim = c(3,n.synapses)) # pre-allocate storage for locations
  
  new.pos <- runif(n = 3,min = -radius,max = radius) #randomly generate candidate position
  synapse.positions[,1] <- new.pos #save spine position
  i = 2 # iterator
  while(TRUE){
    new.pos <- runif(n = 3,min = -radius,max = radius) #randomly generate candidate position
    if(sum(new.pos^2)^0.5 > radius) next; # reject if distance is grader than radius
    if(sum(colSums((synapse.positions[,1:i]-new.pos)^2)^0.5 < 2*spine.radius))next; # reject if new.pos is within 2*spine.radius of another spine
    #from above: colSums((synapse.positions[,1:i]-new.pos)^2)^0.5 is euclidean distance to from current position to each spine 
    #"< 2*spine.radius" returns a logical vector with i elements, the i element is TRUE if placing a synapse at the current position would overlap with any previous synapse
    #taking the sum over all synapses returns 0 if and only if no two synapses overlap
    synapse.positions[,i] <- new.pos #save spine position
    i <- i + 1 #increment counter
    if(i > n.synapses) break; #exit loop when desired number of synapses have been generated
  }
  return(synapse.positions)
}

simulate.spine.input <- function(lambda, active.synapse.fraction, spine.distribution){
  # lambda <- 0.350
  # active.synapse.fraction <- 0.1
  n.synapses <- ncol(spine.distribution)
  active.synapses <- sample(x = 1:n.synapses, replace = FALSE, size = floor(active.synapse.fraction*n.synapses) )
  synapse.is.active <- double(length = n.synapses)
  synapse.is.active[active.synapses] <- 1.0
  spine.inputs <- double(n.synapses)
  for(i in seq(n.synapses)){
    x <- ((colSums((spine.distribution - spine.distribution[,i])^2))^0.5) #distances between i spine to all other spines
    spine.inputs[i] <- sum(synapse.is.active*exp(-x/lambda))
  }
  return(list(n.active.synapses = sum(synapse.is.active), spine.inputs=spine.inputs))
}

mc.spine.inputs <- function(lambda, active.synapse.fraction, spine.distribution, n.sims){
  sim.results <- list(n.active.synapses = NULL, total.spine.input = double(n.sims))
  for(i in 1:n.sims){
    s1 <- simulate.spine.input(lambda = lambda,active.synapse.fraction = active.synapse.fraction,spine.distribution = spine.distribution)
    sim.results$total.spine.input[i] <- sum(s1$spine.inputs)
    sim.results$n.active.synapses <- s1$n.active.synapses      
  }
  return(sim.results)
}