vandrivers <- read.table("vandrivers.txt",header=TRUE)
vandrivers$y <- ts(vandrivers$y,start=1969,frequency=12)
