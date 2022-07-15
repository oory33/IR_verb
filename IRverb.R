## Front Back zero
library(tuneR)
setwd("/Users/ryo/Documents/R/R_study/IR_verb/input")

window <- function(L) {
    0.54 - 0.46 * cos(2 * pi * c(0:(L - 1)) / (L - 1))
} # Hamming
prefixir <- "/Users/ryo/Documents/R/R_study/IR_verb/IR/"
prefix <- "/Users/ryo/Documents/R/R_study/IR_verb/input/"
fnames <- list.files(path = prefix, pattern = "*.wav")
irnames <- list.files(path = prefixir, pattern = "*.wav")

for (filename in fnames) {
    for (irname in irnames) {
        solution <- 0
        solution2 <- 0
        wavl <- readWave(filename)@left
        wavr <- readWave(filename)@right
        srate <- readWave(filename)@samp.rate
        setwd("/Users/ryo/Documents/R/R_study/IR_verb/IR")
        readir <- readWave(irname)@left
        readirr <- readWave(irname)@right
        len <- min(length(readir), length(wavl))
        distance <- len / 100
        wavlch <- c(numeric(len), wavl, numeric(len))
        wavrch <- c(numeric(len), wavr, numeric(len))
        resalt <- numeric(length(wavl))
        resalt2 <- numeric(length(wavr))
        i <- 0
        n <- 0
        while (n <= length(wavlch) - len) {
            if (length(readirr) == 0) {
                readirr <- readir
            } else if (n == 0) {
                da1 <- wavlch[1:len]
                da2 <- wavrch[1:len]
                fda1 <- fft(da1 * window(len), inverse = FALSE)
                fda2 <- fft(da2 * window(len), inverse = FALSE)
                fir <- fft(readir, inverse = FALSE)
                firr <- fft(readirr, inverse = FALSE)
                resalt <- Re(fft(fda1 * fir, inverse = TRUE))
                resalt2 <- Re(fft(fda2 * firr, inverse = TRUE))
            } else {
                das <- wavlch[n:(n + len - 1)]
                dasr <- wavrch[n:(n + len - 1)]
                fdas <- fft(das * window(len), inverse = FALSE)
                fdasr <- fft(dasr * window(len), inverse = FALSE)
                fir <- fft(readir, inverse = FALSE)
                firr <- fft(readirr, inverse = FALSE)
                fress <- Re(fft(fdas * fir, inverse = TRUE))
                fressr <- Re(fft(fdasr * firr, inverse = TRUE))
                ress <- c(numeric(distance * i), fress)
                ressr <- c(numeric(distance * i), fressr)
                resalt <- c(c(resalt, numeric(distance)) + ress)
                resalt2 <- c(c(resalt2, numeric(distance)) + ressr)
            }
            i <- i + 1
            n <- n + distance
        }
        solution <- resalt[(len + 1):(len + length(wavl))]
        solution2 <- resalt2[(len + 1):(len + length(wavr))]
        dat <- normalize(Wave(left = solution, right = solution2, samp.rate = srate, bit = 32, pcm = TRUE), unit = "32", center = TRUE)
        setwd("/Users/ryo/Documents/R/R_study/IR_verb/output")
        writeWave(dat, filename = sprintf("%s_%s.wav", substring(filename, 1, (nchar(filename) - 4)), substring(irname, 1, (nchar(irname) - 4))))
        setwd("/Users/ryo/Documents/R/R_study/IR_verb/input")
    }
}
