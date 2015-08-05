### R code from vignette source 'cacheSweave.Rnw'

###################################################
### code chunk number 1: runSweave (eval = FALSE)
###################################################
## Sweave("foo.Rnw")


###################################################
### code chunk number 2: sleepExample (eval = FALSE)
###################################################
## set.seed(1)
## x <- local({
##     Sys.sleep(10)
##     rnorm(100)
## })
## results <- mean(x)


###################################################
### code chunk number 3: useCacheSweave (eval = FALSE)
###################################################
## library(cacheSweave)
## Sweave("foo.Rnw", driver = cacheSweaveDriver)


###################################################
### code chunk number 4: simpleExpr (eval = FALSE)
###################################################
## x <- 1:100


