library(pepFilter)
library(plyr)

K <- 50  # Total number of 'proteins'
Nk <- rep(100, K)  # Total number of peptides in each protein
N <- sum(Nk)  # Total number of peptides
pi <- 0.25  # Probability of an incorrect assignment

correct.assignment <- rep(1:K, Nk)
U <- runif(N)
tmp1 <- sample(K - 1, size = N, replace = TRUE)
alt.assignment <- tmp1 + as.numeric(tmp1 >= correct.assignment)
actual.assignment <- ifelse(U > pi, correct.assignment, alt.assignment)

X <- sim.byprotein(n = 200, Nk = Nk, s = rep(1, K), tau = rep(1, K), rho = 0.4)

X.list <- dlply(data.frame(index = 1:N, actual.assignment), .(actual.assignment), 
    with, {
        list(data = X[, index], assignment = actual.assignment[1], Conf = rep(100, 
            length(index)))
    })
