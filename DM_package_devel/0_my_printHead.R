

#######################################################
# changing printing options from limma :D
# my.printHead
#######################################################


# install.packages("/home/gosia/R/packages/limma", lib="/home/gosia/R/libraries/3.1.0/")


my.printHead <- function(x)
  #  Print leading 5 elements or rows of atomic object
  #  Gordon Smyth
  #  May 2003.  Last modified 14 April 2009.
{
  if(is.atomic(x)) {
    d <- dim(x)
    if(length(d)<2) which <- "OneD"
    if(length(d)==2) which <- "TwoD"
    if(length(d)>2) which <- "Array"
  } else {
    if(inherits(x,"data.frame")) {
      d <- dim(x)
      which <- "TwoD"
    } else {
      if(is.call(x))
        which <- "Call"
      else {
        if(is.recursive(x))
          which <- "Recursive"
        else
          which <- "Other"
      }
    }
  }
  switch(which,
         OneD={
           n <- length(x)
           if(n > 20) {
             print(x[1:5])
             cat("*",n-5,"more elements ...\n")
           } else
             print(x)
         },
         TwoD={
           n <- d[1]
           if(n > 10) {
             
             x <- x[1:5,]
             
             if(d[2] > 10)
               x <- x[, 1:5]
             
             print(x)
             cat("*",n-5,"more rows ...\n") 
             
             if(d[2] > 10)
             cat("*",d[2] - 5, "more columns ...\n")
             
           } else
             print(x)
         },
         Array={
           n <- d[1]
           if(n > 10) {
             dn <- dimnames(x)
             dim(x) <- c(d[1],prod(d[-1]))
             x <- x[1:5,]
             dim(x) <- c(5,d[-1])
             if(!is.null(dn[[1]])) dn[[1]] <- dn[[1]][1:5]
             dimnames(x) <- dn
             print(x)
             cat("*",n-5,"more rows ...\n")
           } else
             print(x)
         },
         Recursive={
           n <- length(x)
           if(n) {
             i <- names(x)
             if(is.null(i)) i <- seq_len(n)
             if(length(i) >= 2) i <- i[1:2]
             
             
             for (what in i) {
               y <- x[[what]]
               cat("$",what, "\n",sep="")
               
               d <- dim(y)
               if(!is.null(d)){
                 
   
                 if(any(d > 10)) {
                   
                   if(d[1] > 10)
                   y <- y[1:5,]
                   
                   if(d[2] > 10)
                     y <- y[, 1:5]
                   
                   print(y)
                   
                   if(d[1] > 10)
                   cat("*",d[1] - 5,"more rows ...\n") 
                   
                   if(d[2] > 10)
                     cat("*",d[2] - 5, "more columns ...\n")
                   
                 } else
                   print(y)
                 
#                  dimy <- dim(y)
#                  if(any(dimy > 10))
#                    y <- y[1:min(5, dimy[1]), 1:min(5, dimy[2])]
               } else
                  print(y) #Recall(y)
         
             }
             cat("**",n-2,"more elements ...\n")
           }
         },
         Call=,
         Other=print(x)
  )
}



unlockBinding("printHead", as.environment("package:limma"))
assignInNamespace("printHead", my.printHead, ns="limma", envir=as.environment("package:limma"))
assign("printHead", my.printHead, as.environment("package:limma"))
lockBinding("printHead", as.environment("package:limma"))





