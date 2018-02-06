require(tseries, quietly = T)
require(quantmod)

getSymbols("AAPL",src='yahoo',from = as.Date("2003-01-01"), to = as.Date("2008-01-01"))
?getSymbols

r <- read.csv("C:/Users/Umbertojunior/Desktop/data science/Firts Semestr/SDS/hw-3/costituent.txt",head=TRUE,sep=",")

liv<-levels(r$Sector)
new<-data.frame(l,row.names = liv)
?data.frame
l<-sample(CD,3)
CD<-r[r$Sector==liv[1],1]


# mandare solo se con wi-fi
for (i in r$Symbol) {
  getSymbols(i, src='yahoo',from = as.Date("2003-01-01"), to = as.Date("2008-01-01"))
}

dic <- matrix()

dic[,1]<-get.hist.quote(instrument="AAPL", start="2003-01-01", end="2008-01-01",quote= c("Open","Close"), provider="yahoo", drop=TRUE)
plot(dic)
dic[1]
