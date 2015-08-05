hmm_health = initHMM(c('Healthy', 'Fever'),
                     c('normal', 'cold', 'dizzy'),
                     transProbs=matrix(c(0.7,0.3,0.4,0.6),2, byrow=TRUE),
                     emissionProbs=matrix(c(0.5,0.4,0.1,0.1,0.3,0.6),2,byrow=TRUE))

hmm_health
viterbi(hmm_health, c('normal', 'cold', 'dizzy'))