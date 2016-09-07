##### Implementation Begins ####
#This code replicates the FED's paper 

require("BMA")
require("PerformanceAnalytics")
require("ggplot2")
require("quantmod")
require("ROCR")

#Load the data from the FRED

fredVariables<- c('USREC', #NBER Recession variable
                  'DPCERAM1M225NBEA', # Consumption Monthly Change, needs adustments
                  'DSPIC96', # Real disposable personal income
                  'INDPRO', #Industrial Production Index
                  'PERMIT', #Housing Permits
                  'PAYEMS', #NonFarm Payrolls
                  'IC4WSA', #4-week moving average of initial Claims (WEEKLY)
                  'AWHMAN', #Weekly hours, Manufacturing
                  'GS10', #10year Treasury yield
                  'T10Y3MM', #10 Year minus 3 Month. (restricted to 1982)
                  'TEDRATE', #TED Spread (daily)
                  'TWEXOMTH', #TradeWeighted Dollar Index
                  'VIXCLS', #VIX Index, it has to be completed. (Daily)
                  'GS3M', #3 Month Treasury
                  'GS2', #2 Year constant maturity
                  'BAA10YM', #Moodys Corporate Spread BAA (Alternative)
                  'AAA10YM', #Moodys Corporate Spread AAA (Alternative2)
                  'M2SL', #Money Stock
                  'M2REAL', #Money Stock in Real Terms
                  'UNRATE' #Unemployement Rate
)

#Getting all the symbols into the Environment from the FRED.
getSymbols(fredVariables,src='FRED')

### Load the data that is in CSV Format.

### Load the Data
GZSPREAD <- read.csv("~/Upwork/Betasmartz - Black Litterman/FED's Model/Dataset/GZSPREAD.csv")
SPX500 <- read.csv("~/Upwork/Betasmartz - Black Litterman/FED's Model/Dataset/SPX500.csv")
USPMI <- read.csv("~/Upwork/Betasmartz - Black Litterman/FED's Model/Dataset/USPMI.csv")

#Transform the dates into date types
GZSPREAD[,1] <- as.Date(as.character(GZSPREAD[,1]) , format="%d/%m/%Y")
SPX500[,1] <- as.Date(as.character(SPX500[,1]) , format="%d/%m/%Y")
USPMI[,1] <- as.Date(as.character(USPMI[,1]) , format="%d/%m/%Y")

#Transform them into XTS Objects
GZSPREAD <- xts(GZSPREAD[,2],order.by=GZSPREAD[,1])
colnames(GZSPREAD)<- 'GZSPREAD'
SPX500 <- xts(SPX500[,2],order.by=SPX500[,1])
colnames(SPX500)<- 'SPX500'
USPMI <- xts(USPMI[,2],order.by=USPMI[,1])
colnames(USPMI)<- 'USPMI'

## Process the problematic data into monthly uniform data
IC4WSA<-aggregate(IC4WSA, as.yearmon, function(x) tail(x, n=1))
VIXCLS<-aggregate(VIXCLS, as.yearmon, function(x) tail(x, n=1))
TEDRATE<-aggregate(TEDRATE, as.yearmon, function(x) tail(x, n=1))

## Calculate the Varibales that are not observed
levelSlope<-xts(rowMeans(merge(GS3M,GS2,GS10)), order.by=index(GS10))
colnames(levelSlope)<-'levelSlope'
curvature<- xts(2*GS2-(GS3M+GS10),order.by=index(GS3M))
colnames(curvature)<-'curvature'
consumReal <-cumprod(1+DPCERAM1M225NBEA/100)
colnames(consumReal) <- 'consumReal'


## Create the response variables for the analysis
USREC_lag6  <-lag(USREC,-6)
colnames(USREC_lag6) <- 'USREC_lag6'
USREC_lag12 <-lag(USREC,-12)
colnames(USREC_lag12) <- 'USREC_lag12'
USREC_lag24 <-lag(USREC,-24)
colnames(USREC_lag24) <- 'USREC_lag24'

## Recession Matrix
recData <- merge(USREC,USREC_lag6,USREC_lag12,USREC_lag24)

## Create a dataframe of the log 3 month return variables
logReturns3M <- merge(SPX500,M2SL,M2REAL,TWEXOMTH,
                      INDPRO,DSPIC96,consumReal,PERMIT,PAYEMS,IC4WSA,AWHMAN,USPMI)


#Calculate the log of the DF
logReturns3M <- log(logReturns3M)

#Merge the last variable
logReturns3M <- merge(logReturns3M,UNRATE)

## Name the columns
colnames(logReturns3M)<- c('Returns3MSPX500',
                           'Returns3MM2SL',
                           'Returns3MM2REAL',
                           'Returns3MTWEXOMTH',
                           'Returns3MINDPRO',
                           'Returns3MDSPIC96',
                           'Returns3MconsumReal',
                           'Returns3MPERMIT',
                           'Returns3MPAYEMS',
                           'Returns3MIC4WSA',
                           'Returns3MAWHMAN',
                           'Returns3MUSPMI',
                           'Returns3MUNRATE')

#Calculate the 3 Month differences
logReturns3M <- diff(logReturns3M,3)

#Create the DataFrame with the data relevant for the forecast
Data <- merge(DPCERAM1M225NBEA, # Consumption Monthly Change, needs adustments
              DSPIC96, # Real disposable personal income
              INDPRO, #Industrial Production Index
              PERMIT, #Housing Permits
              PAYEMS, #NonFarm Payrolls
              IC4WSA, #4-week moving average of initial Claims (WEEKLY)
              AWHMAN, #Weekly hours, Manufacturing
              GS10, #10year Treasury yield
              T10Y3MM, #10 Year minus 3 Month. (restricted to 1982)
              TEDRATE, #TED Spread (daily)
              TWEXOMTH, #TradeWeighted Dollar Index
              VIXCLS, #VIX Index, it has to be completed. (Daily)
              GS3M, #3 Month Treasury
              GS2, #2 Year constant maturity
              BAA10YM, #Moodys Corporate Spread BAA (Alternative)
              AAA10YM, #Moodys Corporate Spread AAA (Alternative2)
              M2SL, #Money Stock
              M2REAL, #Money Stock in Real Terms
              UNRATE,
              levelSlope,
              curvature,
              consumReal,
              GZSPREAD,
              SPX500,
              USPMI,
              logReturns3M,
              recData )

#Filther the complete cases

Data_clean <- Data[complete.cases(Data),]
working_df<-data.frame(Data_clean)

# #Extract the first dates and the last dates
# initial_dates=which(diff(Data_clean[,"NBER"])==1)-1
# last_dates=which(diff(Data_clean[,"NBER"])==-1)
# up=rep(Inf,length(initial_dates))
# down=rep(-Inf,length(last_dates))
# shade=data.frame(initial_dates,last_dates,up,down)
# 
# aux2=dim(Data_clean)[1]
# # Ploting the data
# ggplot() + 
#         geom_line(aes(x=c(1:aux2), y=SP500NM), color='red',data=fedData_Q_clean)+
#         geom_rect(data=shade, 
#                   mapping=aes(xmin=initial_dates, xmax=last_dates, ymin=up, ymax=down), color='grey', alpha=0.2)


##### Generate the 6 month Forecast Model #####

## Adjust the model and Store it
fed.bic.glm.6M <- bic.glm (USREC_lag6 ~  
                        levelSlope + T10Y3MM + curvature +
                        TEDRATE+ BAA10YM + AAA10YM +
                        Returns3MSPX500+ Returns3MM2SL+ Returns3MM2REAL+
                        Returns3MTWEXOMTH+ Returns3MINDPRO+ Returns3MDSPIC96+ 
                        Returns3MconsumReal+ Returns3MPERMIT+ Returns3MPAYEMS+ 
                        Returns3MIC4WSA+ Returns3MAWHMAN+ Returns3MUSPMI + Returns3MUNRATE +GZSPREAD,
                        data=working_df, glm.family=binomial(link="probit"), OR=500, OR.fix=50,nBest=250)

summary (fed.bic.glm.6M,conditional=T,digits=2)
imageplot.bma(fed.bic.glm.6M,order="probne0")
features6M<- stack(fed.bic.glm.6M$probne0)
features6M[,'name']<-'features6M'
#table(working_df[,"USREC_lag12"],round(predict(fed.bic.glm, newdata = working_df)))

probRec<-predict(fed.bic.glm.6M, newdata = working_df)
working_df[,'probRec6M'] <- probRec

## Ploting the results

#Extract the first dates and the last dates
aux1=1
initial_dates=c(which(diff(Data_clean[,"USREC_lag6"])==1))
last_dates=which(diff(Data_clean[,"USREC_lag6"])==-1)-1
up=rep(Inf,length(initial_dates))
down=rep(-Inf,length(last_dates))
shade=data.frame(initial_dates,last_dates,up,down)

aux2=dim(Data_clean)[1]

# Ploting the data
ggplot() +
        geom_line(aes(x=c(1:aux2), y=probRec), color='red',data=working_df)+
        geom_rect(data=shade,
                  mapping=aes(xmin=initial_dates, xmax=last_dates, ymin=up, ymax=down), color='grey', alpha=0.2) +
        ggtitle("Recession Forecasting 6-Month Ahead") +
        xlab("Time Period") + ylab("Probability of Recession") + labs(fill = "6 Month Ahead Prob")




##### Generate the 12 month Forecast Model #####

## Adjust the model and Store it
fed.bic.glm.12M <- bic.glm (USREC_lag12 ~  
                                   levelSlope + T10Y3MM + curvature +
                                   TEDRATE+ BAA10YM + AAA10YM +
                                   Returns3MSPX500+ Returns3MM2SL+ Returns3MM2REAL+
                                   Returns3MTWEXOMTH+ Returns3MINDPRO+ Returns3MDSPIC96+ 
                                   Returns3MconsumReal+ Returns3MPERMIT+ Returns3MPAYEMS+ 
                                   Returns3MIC4WSA+ Returns3MAWHMAN+ Returns3MUSPMI + Returns3MUNRATE +GZSPREAD,
                           data=working_df, glm.family=binomial(link="probit"), OR=500, OR.fix=50,nBest=250)

summary (fed.bic.glm.12M,conditional=T,digits=2)
imageplot.bma(fed.bic.glm.12M,order="probne0")
features12M<- stack(fed.bic.glm.12M$probne0)
features12M[,'name']<-'features12M'
#table(working_df[,"USREC_lag12"],round(predict(fed.bic.glm, newdata = working_df)))

probRec<-predict(fed.bic.glm.12M, newdata = working_df)
working_df[,'probRec12M'] <- probRec

## Ploting the results

#Extract the first dates and the last dates
aux1=1
initial_dates=c(aux1,which(diff(Data_clean[,"USREC_lag12"])==1))
last_dates=which(diff(Data_clean[,"USREC_lag12"])==-1)-1
up=rep(Inf,length(initial_dates))
down=rep(-Inf,length(last_dates))
shade=data.frame(initial_dates,last_dates,up,down)

aux2=dim(Data_clean)[1]

# Ploting the data
ggplot() +
        geom_line(aes(x=c(1:aux2), y=probRec), color='red',data=working_df)+
        geom_rect(data=shade,
                  mapping=aes(xmin=initial_dates, xmax=last_dates, ymin=up, ymax=down), color='grey', alpha=0.2) +
        ggtitle("Recession Forecasting 12-Month Ahead") +
        xlab("Time Period") + ylab("Probability of Recession") + labs(fill = "12 Month Ahead Prob")


##### Generate the 24 month Forecast Model #####

## Adjust the model and Store it
fed.bic.glm.24M <- bic.glm (USREC_lag24 ~  
                                   levelSlope + T10Y3MM + curvature +
                                   TEDRATE+ BAA10YM + AAA10YM +
                                   Returns3MSPX500+ Returns3MM2SL+ Returns3MM2REAL+
                                   Returns3MTWEXOMTH+ Returns3MINDPRO+ Returns3MDSPIC96+ 
                                   Returns3MconsumReal+ Returns3MPERMIT+ Returns3MPAYEMS+ 
                                   Returns3MIC4WSA+ Returns3MAWHMAN+ Returns3MUSPMI + Returns3MUNRATE +GZSPREAD,
                           data=working_df, glm.family=binomial(link="probit"), OR=500, OR.fix=50,nBest=250)

summary (fed.bic.glm.24M,conditional=T,digits=2)
imageplot.bma(fed.bic.glm.24M,order="probne0")
features24M<- stack(fed.bic.glm.24M$probne0)
features24M[,'name']<- 'features24M'
#table(working_df[,"USREC_lag12"],round(predict(fed.bic.glm, newdata = working_df)))

probRec<-predict(fed.bic.glm.24M, newdata = working_df)
working_df[,'probRec24M'] <- probRec

## Ploting the results

#Extract the first dates and the last dates
aux1=1
initial_dates=c(which(diff(Data_clean[,"USREC_lag24"])==1))
last_dates=which(diff(Data_clean[,"USREC_lag24"])==-1)-1
up=rep(Inf,length(initial_dates))
down=rep(-Inf,length(last_dates))
shade=data.frame(initial_dates,last_dates,up,down)

aux2=dim(Data_clean)[1]

# Ploting the data
ggplot() +
        geom_line(aes(x=c(1:aux2), y=probRec), color='red',data=working_df)+
        geom_rect(data=shade,
                  mapping=aes(xmin=initial_dates, xmax=last_dates, ymin=up, ymax=down), color='grey', alpha=0.2) +
        ggtitle("Recession Forecasting 24-Month Ahead") +
        xlab("Time Period") + ylab("Probability of Recession") + labs(fill = "24 Month Ahead Prob")

## Making a features matrix ##
features_importance=rbind(features6M,features12M,features24M)
order <- c("features6M","features12M","features24M")
features_importance[match(order, features_importance$name),]

## Dodge Bar for feauture importance

# dodged
ggplot(data=features_importance, aes(x=ind, y=values, fill=name)) +
        geom_bar(stat="identity", position=position_dodge())+
        scale_fill_brewer(palette="Paired")+
        theme_minimal() +
        theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
        ggtitle("Summary of Probability of Inclusion" ) +
        xlab("Variable Name") + ylab("Probability of Inclusion")

write.csv(features_importance, "C:/Users/Vitty2/Documents/Upwork/Betasmartz - Black Litterman/FED's Model/featuresummary.csv")

# ggsave("dodged_bowling.jpg")

