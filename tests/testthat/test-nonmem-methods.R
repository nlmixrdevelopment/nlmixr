context("Test Wang2007 data")

ke<-0.5
omega<-0.04
eps<-0.1;

##define function f()
f<- function(eta, time) {
    return(10*exp(-ke*exp(eta)*time))
}

##define the differential function of f()
fp<- function(eta, time) {
    return(10*exp(-ke*exp(eta)*time)*(-ke*time)*exp(eta))
}

test_that("Data are producing the results in the paper", {
    data<-Wang2007_prop_focei;
    data$fp<-fp(data$ETA1, data$TIME)

    data$f<-f(0, data$TIME)

    data$f<-f(data$ETA1, data$TIME)#for FOCE with interaction

    llik<-0

    for (i in 1:10) {
        data1<-data[data$ID==i,]
        ##residual var-cov matrix for additive error
        ## cov<-data1$fp%*%t(data1$fp)*omega+diag(2)*eps
        ##for proportional error
        cov<-data1$fp%*%t(data1$fp)*omega+diag(data1$f**2)*eps
        ginv<-solve(cov)#inverse matrix of residual var-cov matrix
        sec<-t(data1$DV-data1$IPRE+data1$fp*data1$ETA1)%*%ginv%*%
            (data1$DV-data1$IPRE+data1$fp*data1$ETA1)
        frs<-determinant(cov, logarithm=T)$modulus[[1]]
        sum1<-sec+frs
        llik<-llik+sum1
    }
    llik <- round(as.vector(llik), 5)
    expect_equal(llik, 39.20108)

    data<-Wang2007_prop_foce;
    data$fp<-fp(data$ETA1, data$TIME)

    data$f<-f(0, data$TIME)

    ##data$f<-f(data$ETA1, data$TIME)#for FOCE with interaction

    llik<-0

    for (i in 1:10) {
        data1<-data[data$ID==i,]
        ##residual var-cov matrix for additive error
        ## cov<-data1$fp%*%t(data1$fp)*omega+diag(2)*eps
        ##for proportional error
        cov<-data1$fp%*%t(data1$fp)*omega+diag(data1$f**2)*eps
        ginv<-solve(cov)#inverse matrix of residual var-cov matrix
        sec<-t(data1$DV-data1$IPRE+data1$fp*data1$ETA1)%*%ginv%*%
            (data1$DV-data1$IPRE+data1$fp*data1$ETA1)
        frs<-determinant(cov, logarithm=T)$modulus[[1]]
        sum1<-sec+frs
        llik<-llik+sum1
    }
    llik <- round(as.vector(llik), 5)
    expect_equal(llik, 39.20673)

    data<-Wang2007_prop_fo;
    data$fp<-fp(data$ETA1, data$TIME)

    data$f<-f(0, data$TIME)

    ##data$f<-f(data$ETA1, data$TIME)#for FOCE with interaction

    llik<-0

    for (i in 1:10) {
        data1<-data[data$ID==i,]
        ##residual var-cov matrix for additive error
        ## cov<-data1$fp%*%t(data1$fp)*omega+diag(2)*eps
        ##for proportional error
        cov<-data1$fp%*%t(data1$fp)*omega+diag(data1$f**2)*eps
        ginv<-solve(cov)#inverse matrix of residual var-cov matrix
        sec<-t(data1$DV-data1$IPRE+data1$fp*data1$ETA1)%*%ginv%*%
            (data1$DV-data1$IPRE+data1$fp*data1$ETA1)
        frs<-determinant(cov, logarithm=T)$modulus[[1]]
        sum1<-sec+frs
        llik<-llik+sum1
    }
    llik <- round(as.vector(llik), 5)
    expect_equal(llik, 39.21323);

    data<-Wang2007_add_fo;
    data$fp<-fp(data$ETA1, data$TIME)

    data$f<-f(0, data$TIME)

    ##data$f<-f(data$ETA1, data$TIME)#for FOCE with interaction

    llik<-0

    for (i in 1:10) {
        data1<-data[data$ID==i,]
        ##residual var-cov matrix for additive error
        cov<-data1$fp%*%t(data1$fp)*omega+diag(2)*eps
        ##for proportional error
        ## cov<-data1$fp%*%t(data1$fp)*omega+diag(data1$f**2)*eps
        ginv<-solve(cov)#inverse matrix of residual var-cov matrix
        sec<-t(data1$DV-data1$IPRE+data1$fp*data1$ETA1)%*%ginv%*%
            (data1$DV-data1$IPRE+data1$fp*data1$ETA1)
        frs<-determinant(cov, logarithm=T)$modulus[[1]]
        sum1<-sec+frs
        llik<-llik+sum1
    }
    llik <- round(as.vector(llik), 5)
    expect_equal(llik, 0.02580);


    data<-Wang2007_add_foce;
    data$fp<-fp(data$ETA1, data$TIME)

    data$f<-f(0, data$TIME)

    ##data$f<-f(data$ETA1, data$TIME)#for FOCE with interaction

    llik<-0

    for (i in 1:10) {
        data1<-data[data$ID==i,]
        ##residual var-cov matrix for additive error
        cov<-data1$fp%*%t(data1$fp)*omega+diag(2)*eps
        ##for proportional error
        ## cov<-data1$fp%*%t(data1$fp)*omega+diag(data1$f**2)*eps
        ginv<-solve(cov)#inverse matrix of residual var-cov matrix
        sec<-t(data1$DV-data1$IPRE+data1$fp*data1$ETA1)%*%ginv%*%
            (data1$DV-data1$IPRE+data1$fp*data1$ETA1)
        frs<-determinant(cov, logarithm=T)$modulus[[1]]
        sum1<-sec+frs
        llik<-llik+sum1
    }
    llik <- round(as.vector(llik), 5)
    expect_equal(llik, -2.05850);

})


