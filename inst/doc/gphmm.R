## ----options, echo=FALSE, results="hide",message=FALSE, error=FALSE, include=FALSE, autodep=TRUE----
knitr::opts_chunk$set(fig.align="center", cache=FALSE, error=FALSE, message=FALSE, warning=FALSE)
library(gphmm)
library(Biostrings)
library(jsonlite)

## ------------------------------------------------------------------------
paramgphmm = initializeGphmm()
paramgphmm

## ------------------------------------------------------------------------
read = 'ATGCGATGCA'
read

## ------------------------------------------------------------------------
ref = 'ATGTACGATGA'
ref

## ------------------------------------------------------------------------
computegphmm(read = read,
             ref = ref,
             parameters = paramgphmm,
             output = "short")

## ------------------------------------------------------------------------
computegphmm(read = read,
             ref =  ref,
             qv = 20,
             parameters = paramgphmm,
             output = "long")

## ------------------------------------------------------------------------
n = 100

## ------------------------------------------------------------------------
seqs = generateRandomSequences(n = n, meanLen = 100, sdLen = 2,
                               seed = 7373)
writeXStringSet(seqs, 'queries.fasta')
seqs

## ------------------------------------------------------------------------
toCompute = data.frame(query = rep(paste0('s', 1:n), n),
                       ref = rep(paste0('s', 1:n), each = n))
write.table(toCompute, 'toCompute.csv')
head(toCompute)

## ------------------------------------------------------------------------
out = read.table('toCompute_gphmm.csv', stringsAsFactors = F)
head(out)

## ------------------------------------------------------------------------
n = 50

## ------------------------------------------------------------------------
paramgphmm = initializeGphmm()
paramgphmm

## ------------------------------------------------------------------------
seqs = generateRandomSequences(n = n, meanLen = 100, sdLen = 5,
                               prob = paramgphmm$qR, seed = 7373)
seqs

## ------------------------------------------------------------------------
qv = rnorm(n, 20, 5)
qv[qv < 5] = 5
reads = list()
for (i in 1:n){
  reads[[i]] = generateRead(seq = as.character(seqs[i]), 
                            paramgphmm = paramgphmm,
                            qv = qv[i],
                            seed = i)
}
train = c(seqs, DNAStringSet(sapply(reads, '[[', 1)))
names(train) = c(names(train)[1:n], gsub('s', 't', names(train)[1:n]))
csv = data.frame(reads = paste0('t', 1:n), ref = paste0('s', 1:n), qv = qv)

# write files
writeXStringSet(train, 'train.fasta')
write.table(csv, 'train.csv')

## ------------------------------------------------------------------------
plot(density(qv), main = 'Read QV (Phred)')

## ------------------------------------------------------------------------
plot(density(width(seqs)), main = 'Sequence length')

## ------------------------------------------------------------------------
plot(density(width(train[grepl('t', names(train))])), main = 'Read length')

## ------------------------------------------------------------------------
emiTrans = lapply(reads, function(x) computeCounts(x)) 
emiTrans = lapply(lapply(c(1:4), function(i) lapply(emiTrans, '[[', i)), function(x) Reduce('+', x))
names(emiTrans) = c('counts_emissionM', 'counts_emissionD',
                    'counts_emissionI', 'counts_transition')
emiTrans

## ------------------------------------------------------------------------
nucl = c('A', 'C', 'G', 'T')
estimator = fromJSON('train_paramgphmm.json')
names(estimator[['qR']]) = names(estimator[['qX']]) = names(estimator[['qY']]) = 
  colnames(estimator[['pp']]) = rownames(estimator[['pp']]) = nucl
names(estimator[['deltaX']]) = names(estimator[['deltaY']]) = 
  c('intercept', 'slope')
estimator

## ------------------------------------------------------------------------
ll = fromJSON('train_llgphmm.json')
plot(1:length(ll), ll, xlab = 'Iterations', ylab = 'Log Likelihood',
     type = 'l', main = 'Log likelihood')

## ------------------------------------------------------------------------
bias = unlist(mapply('-', estimator, paramgphmm, SIMPLIFY = FALSE))
ppNames = paste0('pp.', paste0(rep(nucl, each = 4), rep(nucl, 4)))
names(bias)[grep('^pp', names(bias))] = ppNames

## ------------------------------------------------------------------------
emission = grepl('A|C|G|T', names(bias))
plot(bias[emission], xaxt = "n", main = 'Emission parameters',
     xlab = '', ylab = 'Bias')
axis(1, at = 1:length(bias[emission]), las = 2, labels = names(bias[emission]))
abline(h = 0)

## ------------------------------------------------------------------------
transition = !emission
plot(bias[transition], xaxt = "n", main = 'Transition parameters',
     xlab = '', ylab = 'Bias')
axis(1, at = 1:length(bias[transition]), las = 2, labels = names(bias[transition]))
abline(h = 0)

## ------------------------------------------------------------------------
# remove generated files
system('rm queries.fasta')
system('rm toCompute_gphmm.csv')
system('rm toCompute.csv')
system('rm train.csv')
system('rm train.fasta')
system('rm train_paramgphmm.json')
system('rm train_llgphmm.json')

