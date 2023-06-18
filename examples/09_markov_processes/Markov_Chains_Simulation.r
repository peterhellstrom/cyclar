library(cyclar)
library(markovchain) # package for Markov chains

# P1 is a matrix with transition frequencies
P1 <- matrix(c(0,7,4,0,0,1,8,0,1,0,0,11,11,0,0,4), 4, 4, byrow = TRUE)
# P is a matrix with transition probabilities
P <- prop.table(P1, 1)
attr(P,"dimnames") <- list(1:4,1:4)
# In this case, row = start position, column = end position
# So matrix have to be transposed...
P <- t(P)
P <- round(P,2)

zms <- markov.sim(P, 101)
phase.tr(zms)

rowSums(phase.tr(zms)$T)
colSums(phase.tr(zms)$T)



statesNames = c("a", "b")
mcA <- new("markovchain",
           transitionMatrix =
             matrix(c(0.7, 0.3, 0.1, 0.9),
                    byrow = TRUE, nrow = 2,
                    dimnames = list(statesNames, statesNames)))

plot(mcA)
summary(mcA)

r.mcA <- rmarkovchain(100, mcA)
outs <- markovchainSequence(n = 100, markovchain = mcA, t0 = "a")
outs2 <- rmarkovchain(n = 20, object = mcA)

markovchainFit(outs)
phase.tr(outs) # compare with my function

# Two-state Markov chain ----

P2 <- markov.2(p=1,q=1)
P2
z2 <- markov.sim(P2,tmax = 100)
markov.plot(z2)
phase.tr(z2)

P2 <- markov.2(p = 0.5, q = 0.5)
P2
z2 <- markov.sim(P2, tmax = 100)
markov.plot(z2)
phase.tr(z2)

# Useful Markov chain-models ----
# model = c("mover-advancer","mover-stayer","random-walk")

P1 <- markov.models(N = 4, q = 0.25, model = "mover-stayer")
P1
colSums(P1)
zms <- markov.sim(P1, tmax = 500, method = "caswell")
markov.plot(zms)

phase.tr(zms)

##
P1 <- markov.models(N = 4, q = 0.25, model = "mover-advancer")
zms <- markov.sim(P1, tmax = 1001, method = "caswell")
markov.plot(zms)
phase.tr(zms)


##
P1 <- markov.models(N = 4, q = 0.25, model = "mover-stayer")
zms <- markov.sim(P1, tmax = 100000, method = "caswell")
phase.tr(zms)
P1

markov.plot(zms, sample = FALSE)
markov.plot(zms)
markov.plot(zms,
            start = 0,
            sample = TRUE,
            samplen = 200)
markov.plot(zms,
            start = 400,
            sample = TRUE,
            samplen = 200)
markov.plot(zms,
            start = 1800,
            sample = TRUE,
            samplen = 500)


# Create sequence of states, Markov chain ----

n <- 5000
z1 <- markov.sim(P, tmax = n, method = "sample")
z2 <- markov.sim(P, tmax = n, method = "caswell")

markov.plot(z1)
markov.plot(z2)

# Create a random, white noise series for comparison
ps <- rowSums(P1) / sum(P1)
#ps <- 1/(1:nrow(P))
zrand <- sample(x = 1:nrow(P), size = n, replace = TRUE, prob = ps)
phase.tr(zrand)$T

Prand <- markov.models(N = 4, q = 0.25, model = "random-walk")
Prand

# Use phases.tr to estimate transition matrix (simulated)
phase.tr(z1)$T
phase.tr(z2)$T
phase.tr(zrand)$T
