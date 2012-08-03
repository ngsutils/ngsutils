########################################################
# Description
#               CP.CI - Clopper-Pearson interval
#
# Usage
#               CP.CI (n, x, CI=0.95)
#   
#   Argument
#               n:          number of observations
#               x:          number of observed successed
#               CI:         confidence interval range
#
#   By Yunlong Liu, January 17, 2011

CP.CI <- function (n, x, n.allele, CI=0.95)
{
        low_ci <- (1-CI)/2;
        high_ci <- (1-low_ci);
        tl <- qbeta(low_ci,x+1,n-x+1);      # calculate the "exact" low boundary
        th <- qbeta(high_ci,x+1,n-x+1);     # calculate the "exact" high boundary

        # adjust the boundaries to the closest nominal level
        res <- max(1/n, 0.5/n.allele); # use a pseudo-count (0.5 as opposed to 1) 
                                       # for allele num to keep resolution a bit tighter.
        tl <- floor(tl/res)*res;
        th <- ceiling(th/res)*res;

        CP.CI <- c(tl, th);
        return (CP.CI);
}