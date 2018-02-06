load("sp500.RData")
require(tseries, quietly = T)
?get.hist.quote

aapl <- get.hist.quote(instrument="AAPL", start="2003-01-01", end="2008-01-01",
                       quote= c("Open","Close"), provider="yahoo", drop=TRUE)

class(aapl)
names(aapl)
head(aapl)

aapl$y <- aapl$Close/aapl$Open

a=1

plot(aapl, main = "Apple Inc. (AAPL)", xlab = "Year")
save(aapl, a,file ='aapl.RData')
load('aapl.RData')

data("EuStockMarkets")
?EuStockMarkets
plot(log(diff(EuStockMarkets)))
plot(diff(log(EuStockMarkets)))
?log
