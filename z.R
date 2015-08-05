pdf("xh.pdf")
hist(rnorm(100))
dev.off()

# Run from shell with
# R CMD BATCH [*.R]