cpdata <- read.csv ("frequency.csv", header = T, sep = ",")
head(ncpdata)
pdf(file="dp-ncp-wm.pdf", width= 7, height=7)
hist(cpdata$MH.,  breaks = 40, col = "steelblue",  freq = T, xlab = "Molecular weight (Da)")
lines(density(cpdata$MH.), col="black",lwd=2)
rug(jitter(cpdata$MH.))
dev.off()



ncpdata <- read.csv ("frequency.csv", header = T, sep = ",")
head(ncpdata)
pdf(file="dp-ncp-len.pdf", width= 7, height=7)
hist(ncpdata$Length,  breaks = 40, col = "#5fb360",  freq = T, xlab = "Peptide length (amino acid)")
lines(density(ncpdata$Length), col="black",lwd=2)
rug(jitter(ncpdata$Length))
dev.off()

cpdata <- read.csv ("frequency.csv", header = T, sep = ",")
head(ncpdata)
pdf(file="dp-ncp-wm.pdf", width= 7, height=7)
hist(cpdata$MH.,  breaks = 40, col = "steelblue",  freq = T, xlab = "Molecular weight (Da)",xlim=c(500,5000))
lines(density(cpdata$MH.), col="black",lwd=2)
rug(jitter(cpdata$MH.))
dev.off()



ncpdata <- read.csv ("frequency.csv", header = T, sep = ",")
head(ncpdata)
pdf(file="dp-ncp-len.pdf", width= 7, height=7)
hist(ncpdata$Length,  breaks = 40, col = "#5fb360",  freq = T, xlab = "Peptide length (amino acid)",xlim=c(5,50))
lines(density(ncpdata$Length), col="black",lwd=2)
rug(jitter(ncpdata$Length))
dev.off()