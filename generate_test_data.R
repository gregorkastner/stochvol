# get data
library(quantmod)
library(xts)

myenv <- new.env()
getSymbols(c("BMW.DE", "MRK.DE", "AAPL", "^GSPC"), from = "2007-01-01", to = "2018-12-31", env = myenv, src = "yahoo")
myenv <- as.list(myenv)
myenv <- lapply(myenv, function (x) x[, 4])
dat <- list(myenv$AAPL["2007-01-01/2007-12-31"],
            myenv$AAPL["2011-05-05/2016-03-25"],
            myenv$GSPC["2011-05-05/2016-03-25"],
            myenv$GSPC["2013-05-01/"],
            myenv$BMW.DE["2009-05-05/2011-03-19"],
            myenv$BMW.DE["2012-08-25/2015-03-18"],
            myenv$MRK.DE["/2009-12-18"],
            myenv$MRK.DE["2012-08-25/2015-03-18"])

library(stochvol)
library(factorstochvol)

data(exrates)

exr <- xts(exrates[, c("AUD", "USD", "MXN")], order.by = as.Date(exrates$date))

dat <- c(dat, list(exr["2004-01-01/2007-12-31"][, "AUD"],
                   exr["2004-01-01/2007-12-31"][, "USD"],
                   exr["2004-01-01/2007-12-31"][, "MXN"],
                   exr["2007-01-25/2010-03-18"][, "USD"],
                   exr["2007-01-25/2010-03-18"][, "MXN"]))
dat <- lapply(dat, na.omit)
saveRDS(dat, file="dat.RDS")

