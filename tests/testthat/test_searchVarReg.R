library(VarReg)
context("checking searchVarReg and seVarReg functions")
data(lidar)
test_that("Checking model 1: no mono, AIC", {
  model1<-searchVarReg(y=lidar$logratio, x=lidar$range,maxknots.m=2, maxknots.v=2,mono.var="none", selection="AIC",maxit=10)
  expect_equal_to_reference(model1, "SearchModel1.rds")
  model1se<-seVarReg(model1$best.model, boot=TRUE, bootreps=2)
  #expect_equal_to_reference(model1se, "Searchse1.rds")
})
test_that("Checking model 2: no mono, BIC", {
  model2<-searchVarReg(y=lidar$logratio, x=lidar$range,maxknots.m=2, maxknots.v=2,mono.var="none", selection="BIC",maxit=10)
  expect_equal_to_reference(model2, "SearchModel2.rds")
  model2se<-seVarReg(model2$best.model, boot=TRUE, bootreps=2)
  #expect_equal_to_reference(model2se, "Searchse2.rds")
})
test_that("Checking model 3: no mono, HQC", {
  model3<-searchVarReg(y=lidar$logratio, x=lidar$range,maxknots.m=2, maxknots.v=2,mono.var="none", selection="HQC",maxit=10)
  expect_equal_to_reference(model3, "SearchModel3.rds")
  model3se<-seVarReg(model3$best.model, boot=TRUE, bootreps=2)
  #expect_equal_to_reference(model3se, "Searchse3.rds")
})
test_that("Checking model 4: inc mono, HQC", {
  model4<-searchVarReg(y=lidar$logratio, x=lidar$range,maxknots.m=2, maxknots.v=2,mono.var="inc", selection="HQC",maxit=10)
  expect_equal_to_reference(model4, "SearchModel4.rds")
  model4se<-seVarReg(model4$best.model, boot=TRUE, bootreps=2)
  #expect_equal_to_reference(model4se, "Searchse4.rds")
})
test_that("Checking model 5: dec mono, HQC", {
  model5<-searchVarReg(y=lidar$logratio, x=lidar$range,maxknots.m=2, maxknots.v=2,mono.var="dec", selection="HQC",maxit=10)
  expect_equal_to_reference(model5, "SearchModel5.rds")
  model5se<-seVarReg(model5$best.model, boot=TRUE, bootreps=2)
  #expect_equal_to_reference(model5se, "Searchse5.rds")
})
