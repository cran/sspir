vandrivers <- read.table("vandrivers.txt",header=TRUE)
vandrivers$y <- stats::ts(vandrivers$y,start=1969,frequency=12)
