# Niche-overlap

# Data from Krebs 1999 Ecological Methodology, Box 13.2, p. 467
sp1 <- c(7,1,286,71,0,0,0)
sp2 <- c(0,0,38,24,30,140,5)
x <- data.frame(sp1,sp2)

# MacArthur-Levins, Pianka modification
# Percentage overlap
# Morisita's measure
# Simplified Morisita index
# Horn's index
# Hurlberts'index

niche.overlap(x)

e.juv <- c(49,13,1,2)
f.juv <- c(67,31,2,9)
g.juv <- c(39,0,3,1)
h.juv <- c(40,15,3,1)
i.juv <- c(42,20,5,3)

e.ad <- c(11,10,0,0)
f.ad <- c(19,13,0,0)
g.ad <- c(36,25,7,3) 
h.ad <- c(29,5,2,4)
i.ad <- c(51,24,7,2)


niche.overlap(data.frame(e.juv, f.juv, g.juv, h.juv, i.juv))
niche.overlap(data.frame(e.ad, f.ad, g.ad, h.ad, i.ad))

z1 <- data.frame(e.juv, f.juv, g.juv, h.juv, i.juv)
z2 <- data.frame(e.ad, f.ad, g.ad, h.ad, i.ad)

z12 <- cbind(rowSums(z1), rowSums(z2))
niche.overlap(z12)

niche.mat(z1)
niche.mat(z2)
