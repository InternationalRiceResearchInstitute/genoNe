#################
# ldNe
#' Function to compute Ne from LD
#' 
#' Its recommended to compute this cromosome-wise
#'  
#'  Parameters:
#'  @param A bp is a vector of distances in basepairs between loci
#'  @param r2 is a vector of r^2 between loci
#'  @param bpcM is the number of kb per cetimorgan, for rice its 250000
#'  
#'  @return a list containing the Ne and its standard error
################

ldNe<- function(bp, r2, bpcM=250000){
  require(doBy)
  dist<- bp/(bpcM*100)#convert distance to morgans
  df<- na.omit(unique(data.frame(r2= r2, dist=dist)))
  df2<- summaryBy(r2~dist, data=df, FUN=mean)
  q<- loess(r2.mean~dist, data=df2)#fit model
  prd<- predict(q)#get expected values of r2
  y<- (1/prd)-1
  x<- q$x*4
  md<- lm(y~0+x)#this fits the model according to Sved 1971
  Ne<- coefficients(md)[1]
  return(list(Ne=Ne, se=sqrt(vcov(md))))
}