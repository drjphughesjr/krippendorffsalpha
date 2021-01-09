
test_that("krippendorffs.alpha(), confint(), and influence() work",
{
	data = matrix(c(1,2,3,3,2,1,4,1,2,NA,NA,NA,
	                1,2,3,3,2,2,4,1,2,5,NA,3,
	                NA,3,3,3,2,3,4,2,2,5,1,NA,
	                1,2,3,3,2,4,4,1,2,5,1,NA), 12, 4)
    set.seed(12)
	fit = krippendorffs.alpha(data, level = "nominal", confint = TRUE, verbose = FALSE,
	                          control = list(bootit = 100, parallel = FALSE))
    alpha.hat = round(fit$alpha.hat, 3)
	names(alpha.hat) = NULL
    expect_equal(alpha.hat, 0.743)
	ci = round(confint(fit), 3)
	names(ci) = NULL
	expect_equal(ci[1], 0.422)
	expect_equal(ci[2], 1)
	ci = round(confint(fit, level = 0.99), 3)
	names(ci) = NULL
	expect_equal(ci[1], 0.336)
	expect_equal(ci[2], 1)
	inf = influence(fit, units = 6)
	dfbeta = round(inf$dfbeta.units, 3)
	names(dfbeta) = NULL
	expect_equal(dfbeta, -0.114)
})

