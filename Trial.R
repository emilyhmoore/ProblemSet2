rm(list=ls())

#############Problem 1#############

#Set the random seed to an arbitrary number
set.seed(100)
#Create a sample satisfying Benford's law
require(BenfordTests) ##Supported by R 3.0.0 or later
x<-rbenf(n=20)

##Takes a vector of numbers and returns the Leemis m statistic. 
##Note that unlike the function that performs this in the BenfordTest
##Package, this function will only calculate for one significant digit. 
benfordtests<-function(x, dist=TRUE, mstat=TRUE, dstat=TRUE){
  require(BenfordTests)##Benford Tests is supported for R 3.0.0 or later only
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
  reslist<-list(rfreq=rfreq, m=m, d=d)
  everything<-all(dist, mstat, dstat)
  if (everything==TRUE) {return(reslist)} 
  if (everything==FALSE & dist==TRUE){print(reslist[1])}
  if (everything==FALSE & mstat==TRUE){print(reslist[2])}
  if(everything==FALSE & dstat==TRUE) {print(reslist[3])}
}

benfordtests(x, dist=FALSE, dstat=FALSE)

require(BenfordTests) ##Note: This is supported by R 3.0.0 or later only

mdist.benftest(x, digits=1) ##Test my own function against the one in the Benford Test package
edist.benftest(x) ##Test my own function against BenfordTests one

##########Print Benfords Function, Problem 2 ############
print.benfords<-function(x){
  require(data.table)
  results<-benfordtests(x) ##Depends on the benfordtests function above
  r<-c(results[[2]], results[[3]]) ##This is a vector with the stats
  ast<-NULL ##empty vector for signif codes
  ##This if else code tells the computer to pick the right number of
  ##asterisks based on significance for the m stat
  asterisks1<-if (results[[2]]>=0.851 & results[[2]]<0.967) {ast="*"
  }else {
    if(results[[2]]>=0.967 & results[[2]]<1.212) {ast="**"} else {
      if(results[[2]]>=1.212) {ast="***"} else {ast=" "}}
  }
  ##This code is the same as above but for the d stat signif codes
  ast2<-NULL
  asterisks2<-if (results[[3]]>=1.212 & results[[3]]<1.330) {ast2="*"
  }else {
    if(results[[3]]>=1.330 & results[[3]]<1.569) {ast2="**"} else {
      if(results[[3]]>=1.569) {ast2="***"} else {ast2=" "}}
  }
  ##This creates a dataframe of the stats and their signifs
  table1<-data.frame(stats=c(r[1], r[2])
                     ,signif=c(ast,ast2)
                     , row.names=c("M", "D")) ##This is a data table of the results
  
  print(table1)##Returns stats and will return asterisks
  print("Significance Codes: <0.01=*** <0.05=** <0.10=* >0.10 wil return blank ") ##Returns code
  
}

print.benfords(x) ##Uses a benford process so should not be sig
y<-seq(1,20)
print.benfords(y) ##Should be sig different from benford process

tester<-function(){
  #Set the random seed to an arbitrary number
  require(BenfordTests) ##Supported for R 3.0.0 or later.
  set.seed(100)
  #Create a sample satisfying Benford's law
  x<-rbenf(n=20)
  #Create a set of data that does NOT satisfy Benford's law
  y<-seq(1, 20)
  
  ##Targetx1 based on relative frequency of significant digits
  ##Easily counted since there were so few numbers
  targetx1<-c(5/20, 6/20, 3/20, 3/20, 1/20, 1/20, 1/20, 0, 0)
  resx1<-as.numeric(benfordtests(x)[[1]])##Results from function
  ##Are they equal?
  truex1<-all(resx1==targetx1) 
  if (truex1==FALSE) {warning("Distribution for Dataset 1 is WRONG!")}
  
  ##Target for the m statistic for x vector
  targetxm<-as.vector(sqrt(20) * abs(6/20- pbenf(1)[2]))
  ##Result from function:
  resx2<-as.vector(benfordtests(x)[[2]])
  ##Are they equal?
  truex2<-resx2==targetxm 
  if (truex2==FALSE){warning("M Stat for Dataset 1 is WRONG!")}
  
  ##Target for the d statistic for x vector
  targetxd<-as.vector(sqrt(20)*sqrt(sum((targetx1-pbenf(1))^2)))
  ##Results for d from function
  resx3<-benfordtests(x)[[3]]
  ##Returns TRUE if the target matches the function results
  truex3<-resx3==targetxd
  if (truex3==FALSE) {warning("D Stat for Dataset 1 is WRONG!")}
  
  ##This is easily calculable by hand because there are so few numbers.
  targety1<-c(11/20, 2/20, rep(1/20, 7))
  resy1<-benfordtests(y)[[1]]
  resy1<-as.numeric(resy1)
  ##Are they equal?
  truey1<-all(resy1==targety1)
  if (truey1==FALSE) {warning("Distribution for Dataset 2 is WRONG!")}

  ##We know from the distribution that 1 occurs 11 times (the most) as the first sig
  ##digit. I subtracted this from the expected value of 1.
  targetym<-as.vector(sqrt(20) * abs(0.55- pbenf(1)[1]))
  ##Result from function:
  resy2<-as.vector(benfordtests(y)[[2]])
  ##Are they equal?
  truey2<-resy2==targetym
  if (truey2==FALSE) {warning("M Stat for Dataset 2 is WRONG!")}
  
  ##Target d for y vector
  targetyd<-as.vector(sqrt(20)*sqrt(sum((targety1-pbenf(1))^2)))
  ##Results for d from function
  resy3<-benfordtests(y)[[3]]
  ##Returns TRUE if the target matches the function results
  truey3<-resy3==targetyd
  if (truey3==FALSE) {warning("D Stat for Dataset 2 is WRONG!")}
  
  ##Vector of Logicals from above
  logvec<-c(truex1, truex2, truex3, truey1, truey2, truey3)
  ##Creating a Matrix of the logicals with column names to return
  ##if there is a problem in the code
  logmat<-as.matrix(truex1)
  logmat<-cbind(logmat, truex2, truex3, truey1, truey2, truey3)
  labelvec<-c("Correct Distribution Dataset 1"
              ,"Correct M Statistic Dataset 1"
              ,"Correct D Statistic Dataset 2"
              ,"Correct Distribution Dataset 2"
              , "Correct M Statistic Dataset 2"
              , "Correct D Statistic Dataset 2")
  colnames(logmat)<-labelvec ##Setting Column Names
  
  ##Creating conditions for returning a passing result.
  testres<-NULL
  if (all(logvec)==TRUE) {testres=TRUE} else {testres=FALSE}
  ##Returns logical indicating whether the test was passed
  print(list(Results=testres))
  ##Returns logical matrix indicating where something failed if it did
  if(testres==FALSE) {print(logmat)}
}

  ##Running the tester function

####NOTE! Tester function will provide warnings when something is wrong.
##For example, if the distribution is calculated improperly, it will have
##errors for both datasets AND the m and d stat as I use the distribution
##to calculate the stats. As such, pay attention to the first warning
##being given first and work down. 

tester()

#############breaking the function 1: Fail to get proper distribution######

benfordtests<-function(x, dist=TRUE, mstat=TRUE, dstat=TRUE){
  require(BenfordTests)##Benford Tests is supported for R 3.0.0 or later only
  fdig<- signifd(x, 1) ##this is the first digits function. Returns first significant digit
  n<-length(fdig) ##Calculate n 
  digfreq<-table(c(fdig, seq(1,9)))-1 ##This is a table including the ##first signif digits in the data and a sequence of 1 to 9. 
  rfreq<-digfreq ##This is the frequency each digit is observed in the actual data
  rfreq_hyp<-pbenf(1) ##pbenf gives the proportion we should get 
  ##in theory if the data matches Benford. This is the same as log (1+1/d)
  m<-sqrt(n)*max(abs(rfreq - rfreq_hyp)) ##As noted rfreq_hyp is the 
  ##same as log (1+1/d), whereas rfreq is the relative frequency of the 
  ##digits for the data in question. The max simply selects whichever is the max
  ##of the digits being seen as per the formula. 
  d<-sqrt(n)*sqrt(sum((rfreq-rfreq_hyp)^2)) ##Same as above but with sum instead for the d stat
  reslist<-list(rfreq=rfreq, m=m, d=d)
  everything<-all(dist, mstat, dstat)
  if (everything==TRUE) {return(reslist)} 
  if (everything==FALSE & dist==TRUE){print(reslist[1])}
  if (everything==FALSE & mstat==TRUE){print(reslist[2])}
  if(everything==FALSE & dstat==TRUE) {print(reslist[3])}
}

tester()
##Tester returns 6 warnings because the distribution was wrong,
##everything was wrong. OOPS! forgot to divide by n when calculating
##rfreq. This will make the distribution be the raw frequency of each
##significant digit, not the relative frequency that we want!

##############Breaking the code 2: Wrong m or d stat#################
benfordtests<-function(x, dist=TRUE, mstat=TRUE, dstat=TRUE){
  require(BenfordTests)##Benford Tests is supported for R 3.0.0 or later only
  fdig<- signifd(x, 1) ##this is the first digits function. Returns first significant digit
  n<-length(fdig) ##Calculate n 
  digfreq<-table(c(fdig, seq(1,9)))-1 ##This is a table including the ##first signif digits in the data and a sequence of 1 to 9. 
  rfreq<-digfreq/n ##This is the frequency each digit is observed in the actual data
  rfreq_hyp<-pbenf(1) ##pbenf gives the proportion we should get 
  ##in theory if the data matches Benford. This is the same as log (1+1/d)
  m<-max(abs(rfreq - rfreq_hyp)) ##As noted rfreq_hyp is the 
  ##same as log (1+1/d), whereas rfreq is the relative frequency of the 
  ##digits for the data in question. The max simply selects whichever is the max
  ##of the digits being seen as per the formula. 
  d<-sqrt(n)*sqrt(sum((rfreq-rfreq_hyp)^2)) ##Same as above but with sum instead for the d stat
  reslist<-list(rfreq=rfreq, m=m, d=d)
  everything<-all(dist, mstat, dstat)
  if (everything==TRUE) {return(reslist)} 
  if (everything==FALSE & dist==TRUE){print(reslist[1])}
  if (everything==FALSE & mstat==TRUE){print(reslist[2])}
  if(everything==FALSE & dstat==TRUE) {print(reslist[3])}
}

tester()
##Throws warning for m stat indicating that something is wrong
##with that calculation. OOPS! Forgot to divide by sqrt(n)



