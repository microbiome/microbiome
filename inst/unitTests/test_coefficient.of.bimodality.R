# For more details on unit testing, see
# http://master.bioconductor.org/developers/how-to/unitTesting-guidelines/

test_coefficient.of.bimodality <- function() {
    #checkEquals(divideBy(4, 2), 2)
    #checkTrue(is.na(divideBy(4, 0)))
    #checkException(expr, msg)
    #obs <- tryCatch(f(), warning=conditionMessage)
    #checkIdentical("unusual condition", obs)
   checkEqualsNumeric(coefficient.of.bimodality(seq(100)),  0.5669112, tolerance=1.0e-3)

}
