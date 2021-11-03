#A pretty full example showing some low level things:

SFHdemo_dust=SFHfunc()
plot(SFHdemo_dust$out[,2:3], log='x', ylim=c(26,17), xlab=BC03lr$Labels$Wavelab,
     ylab='Mag', col=rev(rainbow(26, end=2/3)), pch=16)

SMdemo=SMstarfunc()
SMdemo

SFHdemo_nodust=SFHfunc(tau_birth=0, tau_screen=0)
wave_star=SFHdemo_nodust$flux[,1]

total_atten=sum(c(0,diff(wave_star))*(SFHdemo_nodust$flux[,2]-SFHdemo_dust$flux[,2]))

plot(SFHdemo_nodust$flux, log='xy', ylim=c(1e-20,1e-15), type='l',
     xlab=BC03lr$Labels$Wavelab, ylab='Flux (erg/s/cm^2/Ang)')
lines(SFHdemo_dust$flux, col='grey')
lines(10^seq(5,7,by=0.01), greybody_norm(10^seq(5,7,by=0.01), z=0.1, norm=total_atten),
      col='brown')
Dale_temp=Dale_interp(type='NormSFR')
lines(Dale_Msol$Wave, Dale_temp$Aspec*total_atten, col='red')

wave=sort(c(SFHdemo_dust$flux[,'wave'], Dale_Msol$Wave))
SpecAtten=approxfun(SFHdemo_dust$flux, rule=2)(wave)+
  approxfun(Dale_Msol$Wave, Dale_temp$Aspec, rule=2)(wave/1.1)
lines(wave, SpecAtten, col='darkgreen')

#Some different SFHs with identical stellar mass formed (10^10 Msol):

bursty=SFHfunc(massfunc=massfunc_p6, m1=2, m2=1, m3=0, m4=0, m5=0,
               forcemass=1e10, agemax=1e10)
constant=SFHfunc(massfunc=massfunc_p6, m1=1, m2=1, m3=1, m4=1, m5=1,
                 forcemass=1e10, agemax=1e10)
young=SFHfunc(massfunc=massfunc_p6, m1=0, m2=0, m3=1, m4=1, m5=0, m6=0,
              forcemass=1e10, agemax=1e10)
old=SFHfunc(massfunc=massfunc_p6, m1=0, m2=0, m3=0, m4=0, m5=1, m6=1,
            forcemass=1e10, agemax=1e10)

#SFHs:

plot(bursty$agevec, bursty$SFR, type='l', col='purple', xlim=c(0,1e10), ylim=c(0,5),
     xlab='Age / Yrs', ylab='SFR / Msol/yr')
lines(constant$agevec, constant$SFR, col='blue')
lines(young$agevec, young$SFR, col='darkgreen')
lines(old$agevec, old$SFR, col='red')
legend('top', legend=c('bursty', 'constant', 'young', 'old'), lty=1, 
       col=c('purple', 'blue', 'darkgreen', 'red'))

#Luminosities with default dust attenuation:

wave=bursty$wave_lum
plot(wave, bursty$lum_atten, log='xy', type='l', col='purple', ylim=c(1e-1,1e8), 
     xlab='Wavelength / Ang', ylab='Lsol / Ang')
lines(wave, constant$lum_atten, col='blue')
lines(wave, young$lum_atten, col='darkgreen')
lines(wave, old$lum_atten, col='red')
legend('topright', legend=c('bursty', 'constant', 'young', 'old'), lty=1, 
       col=c('purple', 'blue', 'darkgreen', 'red'))

