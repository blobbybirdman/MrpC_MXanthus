#!/usr/bin/Rscript


peaksNoAnnot <- scan("final-peakDistances-noAnnot.txt", list(0))
peaksNoAnnot <- peaksNoAnnot[[1]]
peaksNoAnnot <- abs(peaksNoAnnot)
cat(median(peaksNoAnnot))
peaksNoAnnotFilt <- peaksNoAnnot[peaksNoAnnot>=-2500 & peaksNoAnnot<=2500]
png("peaksHistNoAnnot.png")
peaksHistNoAnnot = hist(peaksNoAnnotFilt, breaks=25, main="Distribution of distance between replicate peaks", xlab="Distance", ylab="Frequency")
dev.off()

peaksNoFilt <- scan("final-unfilteredPeakDistances-noAnnot.txt", list(0))
peaksNoFilt <- peaksNoFilt[[1]]
peaksUnfiltFilt <- peaksNoFilt[peaksNoFilt>=-250 & peaksNoFilt<=250]
png("peaksHistNoFilt.png")
peaksHistNoFilt = hist(peaksUnfiltFilt, breaks=100, main="Distribution of distance between replicate peaks", xlab="Distance", ylab="Frequency")
dev.off()

peaksUnfiltCumul <- abs(peaksNoFilt)
breaks = seq(0, 150, by=5)
peaksUnfiltCumul.cut = cut(peaksUnfiltCumul, breaks, right=FALSE)
peaksUnfiltCumul.freq = table(peaksUnfiltCumul.cut)
peaksUnfiltCumul.freq = peaksUnfiltCumul.freq/length(peaksUnfiltCumul)
cumfreq = c(0, cumsum(peaksUnfiltCumul.freq)) 
png("peaksUnfiltCumul.png")
plot(breaks, cumfreq, main="Cumulative distribution of distance between replicate peaks", xlab="Distance", ylab="Cumulative Frequency", xlim=c(0, 150), ylim=c(0, .8))
lines(breaks, cumfreq)


#peaksUnfiltCumulHist = hist(peaksUnfiltCumul, main="Cumulative distribution of distance between replicate peaks", xlab="", ylab="")
dev.off()


peaks <- scan("final-peakORFDistances.txt", list(0))
peaks <- peaks[[1]]
peaksFilt <- peaks[peaks>=-2500 & peaks<=2500]
png("peaksHist.png")
peaksHist = hist(peaksFilt, breaks=50, main="Distribution of distance of peaks from nearest TSS", xlab="Distance", ylab="Frequency")
dev.off()

randPeaks <- scan("final-randomORFDistances.txt", list(0))
randPeaks <- randPeaks[[1]]
randPeaksFilt <- randPeaks[randPeaks>=-2500 & randPeaks<=2500]
png("randPeaksHist.png")
randPeaksHist = hist(randPeaksFilt, breaks=50, main="Distribution of distance of randomised peaks from nearest TSS", xlab="Distance", ylab="Frequency")
dev.off()
