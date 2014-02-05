rm(list=ls())

#Set the random seed to an arbitrary number
set.seed(100)
#Create a sample satisfying Benford's law
x<-rbenf(n=20)

##Takes a vector of numbers and returns the Leemis m statistic. 
##Note that unlike the function that performs this in the BenfordTest
##Package, this function will only calculate for one significant digit. 
benfordtests<-function(x){
  fdig<- signifd(x, 1) ##this is the first digits function. Returns first significant digit
  n<-length(fdig) ##Calculate n 
  digfreq<-table(c(fdig, seq(1,9)))-1 ##This is a table including the ##first signif digits in the data and a sequence of 1 to 9. 
  rfreq<-digfreq/n ##This is the frequency each digit is observed in the actual data
  rfreq_hyp<-pbenf(1) ##pbenf gives the proportion we should get 
  ##in theory if the data matches Benford. This is the same as log (1+1/d)
  m<-sqrt(n) * max(abs(rfreq - rfreq_hyp)) ##As noted rfreq_hyp is the 
  ##same as log (1+1/d), whereas rfreq is the relative frequency of the 
  ##digits for the data in question. The max simply selects whichever is the max
  ##of the digits being seen as per the formula. 
  d<-sqrt(n)*sqrt(sum((rfreq-rfreq_hyp)^2)) ##Same as above but with sum instead for the d stat
  return(list(rfreq, m, d))               
}

benfordtests(x)

mdist.benftest(x, digits=1) ##Test my own function against the one in the Benford Test package
edist.benftest(x) ##Test my own function against BenfordTests one

print.benfords<-function(x){
  require(data.table)
  results<-benfordtests(x) ##Depends on the benfordtests function above
  r<-c(results[[2]], results[[3]]) ##This is a vector with the stats
  table1<-data.frame(stats=c(r[1], r[2]), row.names=c("M", "D")) ##This is a data table of the results
  print("Significance Codes: <0.01=*** <0.05=** <0.10=*") ##Returns code
  print(table1)##Returns stats and will return asterisks
}

print.benfords(x) ##Uses a benford process so should not be sig
print.benfords(y) ##Should be sig different from benford process
