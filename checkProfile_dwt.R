require(wavelets)
chipProfile -> read.delim()
chip_dwt -> dwt(chipProfile[,4] + 0.0)
noise_dwt <- dwt(chipProfile[,5] + 0.0) # input

cor_dwt = NA
for (i in 1: length(chip_dwt@W)){cor_dwt[i] <- cor(chip_dwt@W[[i]], noise_dwt@W[[i]]) }