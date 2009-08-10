## zzz.R --- 
## Author          : Claus Dethlefsen
## Created On      : Fri Jan 21 12:34:57 2005
## Last Modified By: Claus Dethlefsen
## Last Modified On: Sun Jan 23 16:04:10 2005
## Update Count    : 7
## Status          : Unknown, Use with caution!
###############################################################################

.First.lib <- function(lib, pkg)
 {
  if((R.version$major == 1) && (as.numeric(R.version$minor) < 9))
    packageDescription <- package.description
  
  cat("\n")
  cat("-------------------------------------------------------------\n")
#  cat(packageDescription("sspir", lib = lib, field="Title"))
#  cat("\n")
  ver  <- packageDescription("sspir", lib = lib, field="Version")
  maint<- packageDescription("sspir", lib = lib, field="Maintainer")
  autho<- packageDescription("sspir", lib = lib, field="Author")
  descr<- packageDescription("sspir", lib = lib, field="Description")
  built<- packageDescription("sspir", lib = lib, field="Built")
  URL  <- packageDescription("sspir", lib = lib, field="URL")
  cat("sspir:",packageDescription("sspir", lib = lib, field="Title"),", version", ver,"\n")
  cat(descr,"\n")
#  cat("sspir",descr,", version", ver,"\n")
#  cat(paste("sspir, version", ver,  "is now loaded\n"))
  cat("Authors:",autho,"\n")
  cat("Maintained by",maint,"\n")
#  cat("Webpage:",URL,"\n")
  cat("\nBuilt:",built,"\n")
  cat("-------------------------------------------------------------\n")

  }

## printline.R --- 
## Author          : Claus Dethlefsen
## Created On      : Fri Jan 21 12:33:33 2005
## Last Modified By: Claus Dethlefsen
## Last Modified On: Fri Jan 21 12:33:34 2005
## Update Count    : 1
## Status          : Unknown, Use with caution!
###############################################################################

"printline" <-
function (s = "-", n = 60) 
    cat(rep(s, n), "\n", sep = "")

## mysolve.R --- 
## Author          : Claus Dethlefsen
## Created On      : Fri Jan 21 12:32:26 2005
## Last Modified By: Claus Dethlefsen
## Last Modified On: Fri Jan 21 12:32:36 2005
## Update Count    : 1
## Status          : Unknown, Use with caution!
###############################################################################

mysolve <- function(M) {
  B <- chol(M)
  B <- solve(B)
  return(B%*%t(B))
}
