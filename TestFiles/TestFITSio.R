require(FITSio)
## Make test image with axis information, write to disk
Z <- matrix(1:15, ncol = 3)
filename <- paste(getwd(), "/test.fits", sep="")
writeFITSim(Z, file = filename, c1 = "Test FITS file",
            crpix = c(1,1), crvaln = c(10, 100), cdeltn = c(8, 2),
            ctypen = c("Distance", "Time"),
            cunitn = c("Furlongs", "Fortnights"))
## Read image back and display
X <-  readFITS(filename)
ax1 <- axVec(1, X$axDat)          # Make axis vector for image
ax2 <- axVec(2, X$axDat)
xlab <- X$axDat$ctype[1]
ylab <- paste(X$axDat$ctype[2], " [", X$axDat$cunit[2], "]", sep = "")
image(ax1, ax2, X$imDat, xlab = xlab, ylab = ylab)
str(X)
X$axDat
X$hdr[1:10]
X$hdr[which(X$hdr=="BITPIX")+1]
### Read back in, modify data, header, and axis data information,
## then write modified version as new file
Z <-  readFITS(filename)
Z$imDat <- Z$imDat + 300
Z$header <- addKwv('SCALE', 1.03, 'test header mod', header=Z$header)
# Z$axDat <- edit(Z$axDat)  # interactive edit
Z$axDat$cdelt[2] <- 20
filename2 <- paste(getwd(), "test.fits", sep="")
writeFITSim(Z$imDat, file=filename2, axDat=Z$axDat, header=Z$header)
## Clean up files to avoid clutter
unlink(filename)
unlink(filename2)
