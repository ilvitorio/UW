#Function Definitions and library requirements

#Libraries required

require("BMA")
require("PerformanceAnalytics")
require("ggplot2")
require("quantmod")
require("ROCR")

#User defined Functions

initialTransform <- function(data){
  #The function makes the initial data transformation from character to numeric
  numVar<-dim(data)
  data[,1]<-as.Date(data[,1],format="%d/%m/%Y")
  data[,-1]<-apply(data[,-1],2,as.numeric)
  return(data)
}

percentileRank <- function(vectorData){
  #The function ranks a generic ranking
  aux<-round(rank(vectorData)/length(vectorData),1)
  return(aux)
}

DFRank <- function(DFxts){
  #The function ranks the data frame and returns a DF of Ranking
  burn_date<-'1999-12'
  initial_index<-paste("/",burn_date,sep="")
  indexer<-paste(burn_date,"/",sep="")
  index_dates<- index(DFxts[indexer])[-1]
  DFRanked<-DFxts
  DFRanked[initial_index,] <- apply(DFxts[initial_index],2,percentileRank)
  
  for (i in index_dates) { 
    #Computes the roll over in the XTS object for taking into account the 
    #date up to the time
    newIndex<-paste("/",as.Date(i),sep="")
    DFtemp<-DFxts[newIndex]
    newDecile<-tail(apply(DFtemp,2,percentileRank),n=1)
    DFRanked[as.Date(i),]<-newDecile
  }
  
  return(DFRanked)
}


trendExtract <- function(tsUni,frequency){
  #Perform a STL decomposition and extrats the trending time series
  STLObj<-stl(ts(tsUni,frequency=frequency),s.window="periodic")
  TObj<-STLObj$time.series[,2]
  return(TObj)
  
}

DFTrend <- function(DFxts,frequency){
  #The function performs a Seasonal-Trend Decomposition on each column 
  #and returns the DF of the trending times series in a ts object
  DFTrendObj <- apply(DFxts,2,trendExtract,frequency=frequency)
  startDate<-as.Date(start(DFxts))
  Year<-as.numeric(format(startDate,"%Y"))
  Month<-as.numeric(format(startDate,"%m"))
  DFTrendObj <- ts(DFTrendObj,start=c(Year,Month),frequency=frequency)
  return(DFTrendObj)
}


tsMMCycle<- function(tsUni,radius){
          #This function returns a dataframe of the coordinates of the maximum and minimum 
          #points in a time series.
          
          #Take all the critical points in the time series
          diffSeries<- diff(sign(diff(tsUni)))
          points_max<-which(diffSeries==-2)+1
          points_min<-which(diffSeries==+2)+1
          
          #Filter the cycle meaningful critical points in the Time Series following a modified
          #Newton search approach
          auxlagged<- cbind(tsUni,
                            lag(tsUni,radius),
                            lag(tsUni,-radius))
          decision<- sign((auxlagged[,2]-auxlagged[,1])*(auxlagged[,1]-auxlagged[,3]))
          #Filtering the critical points into the Meaningful ones.
          points_max_clear<- points_max[decision[points_max+radius] == -1]
          points_min_clear<- points_min[decision[points_min+radius] == -1]
          
          #Construct the Iterator for filling the vector of expansion and contraction. 
          tsCycle <- rep(0,length(tsUni))
          first_point<-min(points_max,points_min)
          
          #Initialization of the Vector
          if (tsUni[1]<tsUni[first_point]) {
            tsCycle[1]<-0
          }else{
            tsCycle[1]<-1
          }
            
          #Fill the time series of 1s and 0s depending on the cycle.
          tsCycle[points_max_clear] <- rep(1,length(points_max_clear))
          tsCycle[points_min_clear] <- rep(-1,length(points_min_clear))
          tsCycle<- cumsum(tsCycle)
          
          return(tsCycle)
}

DFCycles <- function(DFTrend,radius){
            #THe function takes a smoothed time series dataframe and delivers a 1 and 0 for
            #increasing and decreasing periods.
            # 0 For increase
            # 1 For decrease
            # This is consistent with the NBER variable definition
            colnumber<-dim(DFTrend)[2]
            rownumber<-dim(DFTrend)[1]
            DFCycles<-DFTrend
            
            for (i in seq(colnumber) ){
              result<-try(tsMMCycle(DFTrend[,i],radius),silent=TRUE)
              if (class(result)=="try-error") {
                DFCycles[,i]<-rep(NA,rownumber)
              } else{
                DFCycles[,i]<-result
              }
              
            }
            
            startDate<-start(DFTrend)
            frequency<-tsp(DFTrend)[3]
            DFCycles <- ts(DFCycles,start=startDate,frequency=frequency)
            
            return(DFCycles)
  
}

DFStage <- function(tsDec,tsCycle,minT,maxT){
              #The function puts in place the Economic Bussines Cycles for a variable based
              #in its deciles and its cycle state.
  
  #Heritance of the time series variable and Default State
  tsStage<-tsDec
  tsStage<-rep("Equilibrium", length(tsStage))
  
  ### Filtering the Equilibrium to Peak
  index_threshold <- tsDec >= maxT  
  index_composed <- index_threshold  & !(as.logical(as.numeric(tsCycle)))
  tsStage[index_composed]<- c("Equilibrium to Peak")
  
  ### Filtering Peak to Equilibirum
  index_threshold <-  tsDec >= maxT
  index_composed <- index_threshold  & (as.logical(as.numeric(tsCycle)))
  tsStage[index_composed]<- c("Peak to Equilibrium")
  
  ### Filtering Trough to Equilibrium
  index_threshold <-  tsDec <= minT
  index_composed <- index_threshold  & !(as.logical(as.numeric(tsCycle)))
  tsStage[index_composed]<- c("Pit to Equilibrium")
  
  ### Filtering Equilibirum to Trough 
  index_threshold <-  tsDec <= minT
  index_composed <- index_threshold  & (as.logical(as.numeric(tsCycle)))
  tsStage[index_composed]<- c("Equilibrium to Pit")
  ##
  
  return(tsStage)
}


#Begins the coding

## Charge the data and make the initial data transformation
Data_Month_Level_R <- initialTransform(read.csv("~/General Documents/FED Study/FED-master/Datasets/Data_Month_Level_R.csv",stringsAsFactors=FALSE))
Data_Month_Perc_R <- initialTransform(read.csv("~/General Documents/FED Study/FED-master/Datasets/Data_Month_Perc_R.csv",stringsAsFactors=FALSE))
Data_Month_PercCh_R <- initialTransform(read.csv("~/General Documents/FED Study/FED-master/Datasets/Data_Month_PercCh_R.csv",stringsAsFactors=FALSE))
Data_Quarter_Level_R <- initialTransform(read.csv("~/General Documents/FED Study/FED-master/Datasets/Data_Month_PercCh_R.csv",stringsAsFactors=FALSE))

##Chargin the data Second Source
#Data_Month_Level_R <- initialTransform(read.csv("~/Upwork/Betasmartz - Black Litterman/FED's Model/Data_Month_Level_R.csv",stringsAsFactors=FALSE))
#Data_Month_Perc_R <- initialTransform(read.csv("~/Upwork/Betasmartz - Black Litterman/FED's Model/Data_Month_Perc_R.csv",stringsAsFactors=FALSE))
#Data_Month_PercCh_R <- initialTransform(read.csv("~/Upwork/Betasmartz - Black Litterman/FED's Model/Data_Month_PercCh_R.csv",stringsAsFactors=FALSE))
#Data_Quarter_Level_R <- initialTransform(read.csv("~/Upwork/Betasmartz - Black Litterman/FED's Model/Data_Month_PercCh_R.csv",stringsAsFactors=FALSE))



## Transform the variables into XTS DataFrames
Data_Month_Level_ts <- xts(Data_Month_Level_R[,-1] , order.by=Data_Month_Level_R[,1])
Data_Month_Perc_ts <- xts(Data_Month_Perc_R[,-1] , order.by=Data_Month_Perc_R[,1],frequency=12)
Data_Month_PercCh_ts <- xts(Data_Month_PercCh_R[,-1] , order.by=Data_Month_PercCh_R[,1])
Data_Quarter_Level_ts <- xts(Data_Quarter_Level_R[,-1] , order.by=Data_Quarter_Level_R[,1])

## Filter all the cases to work with
Data_Month_Level_clear   <- Data_Month_Level_ts[complete.cases(Data_Month_Level_ts)]
Data_Month_Perc_clear    <- Data_Month_Perc_ts[complete.cases(Data_Month_Perc_ts)]
Data_Month_PercCh_clear  <- Data_Month_PercCh_ts[complete.cases(Data_Month_PercCh_ts)]
Data_Quarter_Level_clear <- Data_Quarter_Level_ts[complete.cases(Data_Quarter_Level_ts)]

## Apply logarithms
Data_Month_Level_log   <- log(Data_Month_Level_clear)
#Data_Quarter_Level_log <- log(Data_Quarter_Level_clear)

## Returns of the data 1-Month, 3-Month ,6-Month and 12-Month from Log
Data_Month_Level_R1 <- diff(Data_Month_Level_log,1)
Data_Month_Level_R3 <- diff(Data_Month_Level_log,3)
Data_Month_Level_R6 <- diff(Data_Month_Level_log,6)
Data_Month_Level_R12 <- diff(Data_Month_Level_log,12)

## Returns of the data 1-Month, 3-Month ,6-Month and 12-Month from Perc
Data_Month_Perc_R1 <- diff(Data_Month_Perc_clear,1)
Data_Month_Perc_R3 <- diff(Data_Month_Perc_clear,3)
Data_Month_Perc_R6 <- diff(Data_Month_Perc_clear,6)
Data_Month_Perc_R12 <- diff(Data_Month_Perc_clear,12)

## Filter again the cases to work with
Data_Month_Level_R1  <-Data_Month_Level_R1[complete.cases(Data_Month_Level_R1)]
Data_Month_Level_R3  <-Data_Month_Level_R3[complete.cases(Data_Month_Level_R3)]
Data_Month_Level_R6  <-Data_Month_Level_R6[complete.cases(Data_Month_Level_R6)]
Data_Month_Level_R12 <-Data_Month_Level_R12[complete.cases(Data_Month_Level_R12)]


Data_Month_Perc_R1  <-Data_Month_Perc_R1[complete.cases(Data_Month_Perc_R1)]
Data_Month_Perc_R3  <-Data_Month_Perc_R3[complete.cases(Data_Month_Perc_R3)]
Data_Month_Perc_R6  <-Data_Month_Perc_R6[complete.cases(Data_Month_Perc_R6)]
Data_Month_Perc_R12  <-Data_Month_Perc_R12[complete.cases(Data_Month_Perc_R12)]



## Rank the data by deciles

#Returns by Deciles
Data_Month_Level_R1_Dec <- DFRank(Data_Month_Level_R1)
Data_Month_Level_R3_Dec <- DFRank(Data_Month_Level_R3)
Data_Month_Level_R6_Dec <- DFRank(Data_Month_Level_R6)
Data_Month_Level_R12_Dec<- DFRank(Data_Month_Level_R12)

Data_Month_Perc_R1_Dec <- DFRank(Data_Month_Perc_R1)
Data_Month_Perc_R3_Dec <- DFRank(Data_Month_Perc_R3)
Data_Month_Perc_R6_Dec <- DFRank(Data_Month_Perc_R6)
Data_Month_Perc_R12_Dec <- DFRank(Data_Month_Perc_R12)

#Percentage Variables by Deciles
Data_Month_Perc_Dec    <- DFRank(Data_Month_Perc_clear)
Data_Month_PercCh_Dec  <- DFRank(Data_Month_PercCh_clear)

#Seasonal Trend Decomposition to find local minima and maxima
Data_Month_Perc_Trend <- DFTrend(Data_Month_Perc_Dec,12)

#Cycle DataFrame
Data_Month_Perc_Cycle <- DFCycles(Data_Month_Perc_Trend,12)


#Stage DataFrame
minimumT<-0.5
maximumT<-0.6
CycleVar <- DFStage(Data_Month_Perc_Trend[,24],Data_Month_Perc_Cycle[,24],minimumT,maximumT)

dataplot<-data.frame(cbind(Data_Month_Perc_Trend[,24]),CycleVar))
dataplot$a<-as.numeric(as.character(dataplot$a))
colnames(dataplot)<-c("a","b")
line<-ggplot(dataplot, aes(x=seq(dim(dataplot)[1]), y=a, color=b)) +ggtitle("Economic Cycle") +xlab("Time")+ylab("Output GDP GAP") +geom_point()




# a<-ts(as.vector(Data_Month_Perc_Dec$ROUTGAP),start=as.Date(start(Data_Month_Perc_Dec$ROUTGAP)),frequency=12)
# stl(a,s.window="periodic")
# 
# radius<-12
# lagged<-cbind(Data_Month_Perc_Trend[,24],lag(Data_Month_Perc_Trend[,24],radius),lag(Data_Month_Perc_Trend[,24],-radius))
# colnames(lagged)<-c("orig","atras","adela")
# testiviri<-cbind(lagged[,2]-lagged[,1],lagged[,1]-lagged[,3],lagged[,1],diff(sign(diff(Data_Month_Perc_Trend[,24]))))
# 
# 
# 
# colnames(testiviri)<- c("var1","var2","orig","detector")