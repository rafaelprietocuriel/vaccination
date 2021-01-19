#### vaccination
#### Code to simulate a SVIR process
#### (Susceptible, Vaccinated, Infected and Recovered)
#### on a network, considering:
#### - the topology of the contact network
#### - The diffusion of anti-vacination views
#### - Some simulated demographics (age and death of individuals)
#### 
#### - It produces:
#### - The size of the anti-vaccine community
#### - The time for anti-vaccine views to propagate
#### - The size of the recovered individuals
#### - The time for the virus to stop spreading
#### - The years of life lost
#### - The number of casualties

require(igraph)
N <- 5000 #### number of nodes

#### parameters
p <- 0.05 ### probability of infecting a person each day
TotSims <- 50 ### number of simulations

#### functions
#### create a proximity network
GeoN <- function(N, delta = 2){
  x <- runif(N); y <- runif(N)
  return(graph_from_adjacency_matrix((as.matrix(dist(cbind(x, y), 
                                                     upper = T, diag = T)) < delta)*1, 
                                     mode = "undirected",
                                     diag = F))
  
}

#### dropping random edges
#### a function with input: Adj - the adjacency matrix
####                        rho - the % to keep
dropping <- function(Adj, rho = 1){
  #Nearby <- ((Adj %*% Adj)>0)*1
  Adj[lower.tri(Adj)] <- 0
  n <- dim(Adj)[1]
  S <- (matrix(runif(n^2), ncol = n) * Adj) > rho
  NewAdj <- Adj
  NewAdj[S] <- 0 ### remove
  NewAdj <- ((NewAdj + t(NewAdj))>0)*1
  diag(NewAdj) <- 0
  return(NewAdj)
}

#### rewiring random
#### a function with input: Adj - the adjacency matrix
####                        rho - the % to rewire
rewiring <- function(Adj, rho = 1){
  Adj[lower.tri(Adj)] <- 0
  n <- dim(Adj)[1]
  S <- (matrix(runif(n^2), ncol = n) * Adj) > rho
  NewAdj <- Adj
  NewAdj[S] <- 0 ### remove
  S <- which(S, arr.ind = T)
  for(k in 1:dim(S)[1]){ ### rewire
    u <- 1+(runif(1) < 0.5)*1 ### rewire first or second entry
    v <- S[k, 3-u]
    w <- sample.int(n, 1)
    NewAdj[v,w] <- 1; NewAdj[w,v] <- 1
  }
  NewAdj <- ((NewAdj + t(NewAdj))>0)*1
  diag(NewAdj) <- 0
  return(NewAdj)
}

#### create a population with age and death
{
  Age <- 100*rbeta(N, 2, 3)
  Death <-  100*rbeta(N, 5, 2)
  ind <- Age > Death
  while(sum(ind) > 0){
    Age[ind] <- 100*rbeta(sum(ind), 2, 3)
    Death[ind] <- 100*rbeta(sum(ind), 5, 2)
    ind <- Age > Death
  }
}  

#### emtpy Data Frame for simulation results
Strategy <- data.frame(vaccineRate = c(),
                       persuasiveness = c(),
                       antivTime = c(),
                       TotalAntiV = c(),
                       LostLife = c(),
                       recovered = c(),
                       casualties = c(),
                       time = c())

#### create adjacency matrix for social infection
#### also return centrality and node degree
#### Topology of the networks
G <- barabasi.game(N, 
                   power = 1,
                   m = 3, 
                   out.dist = NULL, 
                   out.seq = NULL, 
                   out.pref = FALSE, 
                   directed=FALSE)
G <- sample_smallworld(1,
                       size=N,
                       nei=3,p=0.1)
G <- GeoN(N, 0.019)

#### create adjacency and degree matrices
Adj <- as.matrix(get.adjacency(G))
Deg <- degree(G)
Bet <- betweenness(G, directed = F)
RewiredAdj <- Adj
RewiredDeg <- Deg


#### rewire and dropping edges
#### rewire 30% of the edges
NewAdj <- rewiring(Adj, 0.7)
G <- graph_from_adjacency_matrix(NewAdj, 
                                   mode = "undirected")
Adj0.7 <- as.matrix(get.adjacency(G))
Deg0.7 <- degree(G)
Bet0.7 <- betweenness(G, directed = F)

#### Dropping 40% of the edges  
NewAdj <- dropping(Adj, 0.6)
G <- graph_from_adjacency_matrix(NewAdj, 
                                   mode = "undirected")
Adj0.4 <- as.matrix(get.adjacency(G))
Deg0.4 <- degree(G)
Bet0.4 <- betweenness(G, directed = F)
  

#### Contact network 
#### it can be either the dropping edges, the rewired or the same
RewiredAdj <- Adj0.7
RewiredDeg <- Deg0.7


#### fixed vaccination rate or fixed antivaccination persuasiveness
for(sims in 1:TotSims){
  vaccineRate <- 0.3
  # vaccineRate <- runif(1)
  antiv <- runif(1)
  # antiv <- 0.2
  #### first the antivaccination part
  {
      A <- ((1:N) %in% sample.int(N, N/80))*1
      E <- (1 - A) * apply(cbind(Adj[, A >= 1], rep(0, N)), 1, max)
      S <- 1 - A - E
      D <- 0*E ### decided after being exposed to antivax
      
      ### Y and N and just cosmetic, Y said Y to AV
      ### N said no to AV 
      ### timestamp
      tantiv <- 0 
      
      #### time vars
      At <- sum(A>0)
      St <- sum(S>0)
      Et <- sum(E>0)
      Dt <- sum(D>0)
      Nt <- 0
      Yt <- 0
      
      nExposed <- sum(E)
      while(nExposed > 0){#### there are exposed people
        NewA <- runif(N) * E > (1-antiv)
        A[NewA] <- 1
        D[which(E - NewA>0)] <- 1
        E <- (1 - A - D) * apply(cbind(Adj[, A >= 1], rep(0, N)), 1, max)
        S <- 1 - A - E - D
        
        #### time
        tantiv <- tantiv + 1
        
        #### update times
        At <- c(At, sum(A>0))
        Dt <- c(Dt, sum(D>0))
        Et <- c(Et, sum(E>0))
        St <- c(St, sum(S>0))
        Yt <- c(Yt, sum(NewA))
        Nt <- c(Nt, nExposed - sum(NewA))
        nExposed <- sum(E)
      }
    }
  #SIR dynamics removing antivacc
  {
    I <- ((1:N) %in% sample.int(N, N/80))*1
      
    ### vaccinate, removing infected and antivaccers
    ### vaccinate most central nodes
    V <- ((1:N) %in% order(RewiredDeg, decreasing = T)[1:(N*vaccineRate)])*(1-I)*(1-A)
    
    #### vaccinate old people
    V <- ((1:N) %in% order(Age, decreasing = T)[1:(N*vaccineRate)])*(1-I)*(1-A)
    
    #### vaccinate young people
    V <- ((1:N) %in% order(Age, decreasing = F)[1:(N*vaccineRate)])*(1-I)*(1-A)
    
    #### vaccinate random nodes
    V <- ((1:N) %in% sample.int(N)[1:(N*vaccineRate)])*(1-I)*(1-A)
    S <- 1 - I - V
      
    ### recoveries
    R <- 0 * I
      
    ### timestamp
    t <- 0 
      
    #### time vars
      It <- sum(I>0)
      St <- sum(S>0)
      Rt <- sum(R>0)
      Vt <- sum(V>0)
      while(sum(I) > 0){#### there are infected people
        I[I>0] <- I[I>0] + 1
        
        #### infect
        #### they have an infected person and are succeptible
        NewI <- runif(N) * S * apply(cbind(RewiredAdj[, I >= 1], rep(0, N)), 1, max) > (1-p)
        I[NewI] <- 1
        S[NewI] <- 0
        
        #### I to R
        R[I > 14] <- 1
        I[I > 14] <- 0
        
        #### time
        t <- t + 1
        
        #### update times
        It <- c(It, sum(I>0))
        St <- c(St, sum(S>0))
        Rt <- c(Rt, sum(R>0))
        vt <- c(Vt, sum(Vt>0))
      }
      #### Deaths
      #### computed only for recovered
      Deaths <- runif(N)*200 < Age * R
      LostYears <- sum((Death - Age)[Deaths])/N
      cat("\n Antivaccine Rate: ", antiv, "  AntiVax ", sum(A), " Av Steps ", tantiv, "\n")
      cat("Vaccine Rate: ", vaccineRate, "  Lost years ",LostYears, " SIR time ",t, "\n")
      cat(" - - - - - - - - - - - - - - - - - - - - \n")
      Strategy <- rbind(Strategy, 
                        data.frame(vaccineRate = vaccineRate,
                                   persuasiveness = antiv,
                                   antivTime = tantiv,
                                   TotalAntiV = sum(A)/N,
                                   LostLife = LostYears,
                                   recovered = sum(R)/N,
                                   casualties = Deaths/N,
                                   time = t))
    }
  }
plot(Strategy)
