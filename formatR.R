### R code from vignette source 'formatR.Rnw'

###################################################
### code chunk number 1: no-comment
###################################################
if (require('knitr')) opts_chunk$set(comment=NA)


###################################################
### code chunk number 2: usage-tidy-source
###################################################
library(formatR)
usage(tidy.source, width = .73)


###################################################
### code chunk number 3: setup (eval = FALSE)
###################################################
## options(keep.blank.line = FALSE)


###################################################
### code chunk number 4: width-spec
###################################################
options(width = 80)


###################################################
### code chunk number 5: tidy-source-ex
###################################################
library(formatR)
## use the 'text' argument
src = c("    ## comments are retained;",
    "# a comment block will be reflowed if it contains long comments;", 
    "#' roxygen comments will not be wrapped in any case",
    "1+1", "if(TRUE){", 
    "x=1  # inline comments", "}else{", "x=2;print('Oh no... ask the right bracket to go away!')}", 
    "1*3 # one space before this comment will become two!", "2+2+2    # 'short comments'", "   ", "lm(y~x1+x2)  ### only 'single quotes' are allowed in comments", 
    "'a character string with \t in it'",
    "1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1  ## comments after a long line", paste("#' here is a", 
        paste(rep("long", 20), collapse = " "), "roxygen comment"))


###################################################
### code chunk number 6: source-code
###################################################
cat(src, sep = "\n")


###################################################
### code chunk number 7: tidy-replace-assign
###################################################
tidy.source(text = src[1:8], replace.assign=TRUE)


###################################################
### code chunk number 8: blank-lines
###################################################
## note the 11th line [an empty line] was discarded
tidy.source(text = src, keep.blank.line = FALSE)


###################################################
### code chunk number 9: reindent-spaces
###################################################
tidy.source(text = src, reindent.spaces = 2)


###################################################
### code chunk number 10: left-brace-newline
###################################################
tidy.source(text = src[1:10], left.brace.newline = TRUE)


###################################################
### code chunk number 11: discard-comments
###################################################
tidy.source(text = src, keep.comment = FALSE)


###################################################
### code chunk number 12: lyx-doc (eval = FALSE)
###################################################
## system.file('doc', 'formatR.Rnw', package='formatR')


###################################################
### code chunk number 13: write-bib
###################################################
if (require('knitr'))
write_bib(c('Rd2roxygen', 'knitr', 'formatR'), file = 'formatR.bib')


