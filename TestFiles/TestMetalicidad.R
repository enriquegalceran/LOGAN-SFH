# Comprobando metalicidades

setwd("~/Documents/GitHub/LOGAN-SFH")
source("FuncionesAuxiliaresYLibrerias.R")

agevec = seq(0, 13.8e9, 1e8)
massvec = massfunc_snorm_burst(agevec, mpeak=12, mburst=5, mburstage=4)
Zp2 = Zfunc_p2(agevec, Z1=0.03, Z2=0.01)
Zlin = Zfunc_massmap_lin(agevec, Zstart=0.01, Zfinal=0.03, massfunc=massfunc_snorm_burst, mpeak=12, mburst=5, mburstage=4)
Zbox = Zfunc_massmap_box(agevec, Zstart=0.01, Zfinal=0.03, massfunc=massfunc_snorm_burst, mpeak=12, mburst=5, mburstage=4)

twoord.plot(agevec, Zp2, agevec, massvec, type="l", lcol="blue", rcol="black",
            main="massfunc_snorm_burst")
lines(agevec, Zlin, col="red")
lines(agevec, Zbox, col="green")
legend("bottomleft", legend=c("Zp2", "Zli", "Zbox"), lty=c(1, 1, 1), col=c("blue", "red", "green"))
legend("bottom", legend="massvec", col="black", lty=1)


Stars = SFHfunc(
  z = 0,
  stellpop = "EMILES",
  massfunc = massfunc_snorm_burst, mpeak=12, mburst=5, mburstage=4,
  emission = FALSE,
  Z=Zfunc_massmap_box,
  Zstar=0.01, Zfinal=0.03
)

plot(Stars$agevec, Stars$SFR, type="l")

plot(Stars$agevec, Stars$massvec, type="l")
plot(Stars$agevec, Stars$Zvec, type="l")
data("EMILES")
plot(EMILES$Age, EMILES$AgeWeights, type="l")

