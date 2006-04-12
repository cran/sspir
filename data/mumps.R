mumps <- read.table("mumps.txt",header=TRUE)
mumps <- ts(mumps,start=1928,end=c(1972,6),frequency=12)
