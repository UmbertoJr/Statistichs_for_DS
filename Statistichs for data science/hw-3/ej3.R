require(tseries, quietly = T)


IT = c("ACN", "ATVI", "ADBE", "AKAM", "ADS", "ADI", "AAPL", "AMAT", "ADSK", "ADP", "CA", "CSCO", "CTXS", "CTSH", "EBAY", "EA", "FFIV", "FIS", "FISV", "GPN", "HRS", "HPQ", "INTC", "IBM", "INTU", "JNPR", "KLAC", "LRCX", "LLTC", "MCHP", "MU", "MSFT", "MSI", "NTAP", "NFLX", "NVDA", "ORCL", "PAYX", "QCOM", "RHT", "STX", "SWKS", "SYMC", "TXN", "TSS", "VRSN", "WDC", "XRX", "XLNX", "YHOO")
HE = c("ABT", "AET", "A", "AGN", "ALXN", "ABC", "AMGN", "ANTM", "BCR", "BAX", "BDX", "BIIB", "BSX", "BMY", "CAH", "HSIC", "CELG", "CNC", "CERN", "CI", "DVA", "XRAY", "EW", "ENDP", "ESRX", "GILD", "HOLX", "HUM", "ILMN", "ISRG", "JNJ", "LH", "LLY", "MCK", "MDT", "MRK", "MTD", "MYL", "PDCO", "PKI", "PRGO", "PFE", "DGX", "REGN", "STJ", "SYK", "COO", "TMO", "UNH", "UHS", "VAR", "VRTX", "WAT", "ZBH")
MA = c("APD", "ALB", "AMAT", "AVY", "BLL", "DOW", "DD", "EMN", "ECL", "FMC", "FCX", "GWW", "IP", "IFF", "MLM", "MON", "MOS", "NEM", "NUE", "PPG", "PX", "SEE", "SHW", "VMC")
EN = c("AES", "LNT", "APC", "APA", "BHI", "COG", "CNP", "CHK", "CVX", "XEC", "CMS",  "COP", "DVN", "DTE", "DUK", "EOG", "EQT", "ES", "XOM", "FE", "FTI", "HAL", "HP", "HES", "MRO", "MUR", "NOV", "NFX", "NBL", "OXY", "OKE", "PXD", "RRC", "SLB", "SRE", "SWN", "TSO", "RIG", "VLO", "WEC", "WMB", "XEL")
stocks = c(IT, HE, MA, EN)

n = length(stocks)
alpha = 0.05
eps = 0.45

data = matrix(nrow = 1257, ncol = length(stocks))
colnames(data) = stocks

#https://en.wikipedia.org/wiki/List_of_S%26P_500_companies

for (i in 1:length(stocks)){
  s = get.hist.quote(instrument = stocks[i], start = "2003-01-01", end = "2008-01-01", quote = c("Open", "Close"), provider = "yahoo", drop = T)
  data[,i] = diff(log(s$Close))
  print(stocks[i])
  }


# Transformation
trans = function(x) 1/2*log((1+x)/(1-x))
trans.inv = function(t) (exp(2*t)- 1)/(exp(2*t) + 1)
  

# Let's see the correlations:

R = cor(data, method = "pearson")
Z = trans(R)

#Bootstrap
B = 4000
delta = rep(NA, B)

for (i in 1:B){
  data.boot = data[sample(nrow(data), replace = T),]
  R.boot = cor(data.boot)
  delta[i] = sqrt(n)*max(abs(R.boot-R))
}


# The ecdf is
F.hat = ecdf(delta)
plot(F.hat)
t.alpha = quantile(F.hat, alpha)
 
# The confidence interval limits are
L = R-t.alpha/sqrt(n)
U = R+t.alpha/sqrt(n)

# Now we have to put edges when the following condition holds

put.edges = (L > eps)*1
diag(put.edges) = 0
put.edges[put.edges>0]

library(igraph)

# Create the graph
g = graph.adjacency(put.edges) 
for(i in 1:length(stocks)){
  if (V(g)$name[i] %in% IT  == T ) V(g)[i]$color = "red"
  if (V(g)$name[i] %in% HE  == T ) V(g)[i]$color = "pink"
  if (V(g)$name[i] %in% MA  == T ) V(g)[i]$color = "gray"
  if (V(g)$name[i] %in% EN  == T ) V(g)[i]$color = "green"
}

plot(g, edge.color="black", vertex.size=7, layout = layout.fruchterman.reingold)
