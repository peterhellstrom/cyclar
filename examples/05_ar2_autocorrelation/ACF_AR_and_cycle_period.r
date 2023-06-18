library(cyclar)

z <- ar2.acf.test()
z

ar2.acf.test(k=3, npoints=50, n=100, nrepl=100, k.per=c("3","6","7","9"))

ar2.acf.test(k=3.5, npoints=50, n=100, nrepl=100, k.per=c("3","6","7","14"))

ar2.acf.test(k=4, npoints=50, n=100, nrepl=100, k.per=c("3","4","5","8","9"))

ar2.acf.test(k=4.5, npoints=50, n=100, nrepl=200, k.per=c("3","4","5","8","9"))
