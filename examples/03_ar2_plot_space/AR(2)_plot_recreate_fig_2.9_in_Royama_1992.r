library(cyclar)

# Fig. 2.9, p. 63 in Royama (1992)

x <- c(0.4, 0.6, -0.8, -0.8, -1, -1, 1.5, 1.5, -1.3, -2.5, 1.5, 1, 0.47) # (1+a1)
y <- c(0.3, 0.6, 0.1, 0.6, -0.8, -1.05, -0.9, -1.05, -1.01, -0.5, -0.75, 0, 0.47) # a2

ar2.plot(k = NULL, v = NULL, xlim = c(-2.5, 2.5), ylim = c(-1.5, 1.5), main = "Fig 2.9 in Royama (1992), p. 63", text.lab = FALSE, par.name = "alpha")
points(x, y, col = 2, pch = 16)
text(x + .15, y, paste("(", 1:13, ")", sep = ""), cex = 0.7)

ar2.period(a1 = x - 1, a2 = y)
