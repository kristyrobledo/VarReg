# building varreg.multi
mean.vars<-c(1,2)
var.vars<-c(1,2)
mean.model<-c('linear', 'linear')
var.model<-c('linear', 'linear')
knots.m<-c(2,2)
knots.v<-NULL
degree=2
var.model<-c('linear', 'semi')
knots.v<-c(NA,2)
mean.model<-c('linear', 'semi')
knots.m<-c(NA, 2)

mean.model="zero" #mean.vars = NULL,
var.model=c("linear", "linear", "linear")
var.vars = c(1,2,3)
knots.v = NULL


########################################################################
## testing varreg.multi

#library(VarReg)
library(palmerpenguins)
library(tidyverse)

cc<-penguins %>%
  na.omit()
y<-cc$body_mass_g
x<-as.data.frame(cbind(cc$bill_length_mm, cc$flipper_length_mm,cc$bill_depth_mm,
                       cc$bill_length_mm*sqrt(2), cc$flipper_length_mm*sqrt(2),cc$bill_depth_mm*sqrt(2)))
colnames(x) <-c("bill length mm", "flipper length mm","bill depth mm",
                "bill_length_mm*sqrt(2)", "flipper_length_mm*sqrt(2)","bill_depth_mm*sqrt(2)")
colnames(x)

##################################################################################
# zero mean, constant variance
m1<-semiVarReg.multi(cc$body_mass_g, x,
                           mean.model="zero", #mean.vars = NULL,
                           var.model=c("constant"), #var.vars = NULL,
                           maxit=1000)
m1[-21]
# zero mean, one variance (linear)
m2<-semiVarReg.multi(cc$body_mass_g, x,
                     mean.model="zero", #mean.vars = NULL,
                     var.model=c("linear"), var.vars = c(1),
                     maxit=1000)
m2[-21]

# zero mean, one variance (semi)
m3<-semiVarReg.multi(cc$body_mass_g, x,
                     mean.model="zero", #mean.vars = NULL,
                     var.model=c("semi"), var.vars = c(1), knots.v = 3,
                     maxit=1000)
m3[-21]

# zero mean, 3 variance (linear)
m4<-semiVarReg.multi(cc$body_mass_g, x,
                     mean.model="zero", #mean.vars = NULL,
                     var.model=c("linear", "linear", "linear"), var.vars = c(1,2,3),knots.v = NULL,
                     maxit=10)
m4[-21]
# zero mean, 3 variance (semi 2 linear 1)
m5<-semiVarReg.multi(cc$body_mass_g, x,
                     mean.model="zero", #mean.vars = NULL,
                     var.model=c("semi", "semi", "linear"), var.vars = c(1,2,3),knots.v = c(2,3),
                     maxit=10)
m5[-21]

##################################################################################
# constant mean constant variance
m6<-semiVarReg.multi(cc$body_mass_g, x,
                     mean.model="constant", #mean.vars = NULL,
                     var.model=c("constant"), #var.vars = NULL,
                     maxit=1000)
m6[-21]

# constant mean one variance (linear)
m7<-semiVarReg.multi(cc$body_mass_g, x,
                     mean.model="constant", #mean.vars = NULL,
                     var.model=c("linear"), var.vars = c(1),
                     maxit=1000)
m7[-21]


# constant mean one variance (semi)
m8<-semiVarReg.multi(cc$body_mass_g, x,
                     mean.model="constant", #mean.vars = NULL,
                     var.model=c("semi"), var.vars = c(1), knots.v = 3,
                     maxit=1000)
m8[-21]

# constant mean 3 variance (linear)
m9<-semiVarReg.multi(cc$body_mass_g, x,
                     mean.model="constant", #mean.vars = NULL,
                     var.model=c("linear", "linear", "linear"), var.vars = c(1,2,3),knots.v = NULL,
                     maxit=10)
m9[-21]
# constant mean 3 variance (semi 2 linear 1)

m10<-semiVarReg.multi(cc$body_mass_g, x,
                     mean.model="constant", #mean.vars = NULL,
                     var.model=c("semi", "semi", "linear"), var.vars = c(1,2,3),knots.v = c(2,3),
                     maxit=10)
m10[-21]

############################################################################
# linear mean, constant variance
m11<-semiVarReg.multi(cc$body_mass_g, x,
                     mean.model="linear", mean.vars = 2,
                     var.model=c("constant"), #var.vars = NULL,
                     maxit=1000)
m11[-21]

# linear mean, one variance (linear)
m12<-semiVarReg.multi(cc$body_mass_g, x,
                      mean.model="linear", mean.vars = 2,
                     var.model=c("linear"), var.vars = c(1),
                     maxit=1000)
m12[-21]

# linear mean, one variance (semi)
m13<-semiVarReg.multi(cc$body_mass_g, x,
                      mean.model="linear", mean.vars = 2,
                     var.model=c("semi"), var.vars = c(1), knots.v = 3,
                     maxit=1000)
m13[-21]
# linear mean, 3 variance (linear)
m14<-semiVarReg.multi(cc$body_mass_g, x,
                      mean.model="linear", mean.vars = 2,
                     var.model=c("linear", "linear", "linear"), var.vars = c(1,2,3),knots.v = NULL,
                     maxit=10)
m14[-21]

# linear mean, 3 variance (semi 2 linear 1)
m15<-semiVarReg.multi(cc$body_mass_g, x,
                      mean.model="linear", mean.vars = 2,
                      var.model=c("semi", "semi", "linear"), var.vars = c(1,2,3),knots.v = c(2,3),
                      maxit=10)
m15[-21]

###########################################################################
# semi mean, constant variance
m16<-semiVarReg.multi(cc$body_mass_g, x,
                      mean.model="semi", mean.vars = 2, knots.m=3,
                      var.model=c("constant"), #var.vars = NULL,
                      maxit=10)
m16[-21]
# semi mean, one variance (linear)
m17<-semiVarReg.multi(cc$body_mass_g, x,
                      mean.model="semi", mean.vars = 2, knots.m=3,
                      var.model=c("linear"), var.vars = c(1),
                      maxit=10)
m17[-21]

# semi mean, one variance (semi)
m18<-semiVarReg.multi(cc$body_mass_g, x,
                      mean.model="semi", mean.vars = 2, knots.m=3,
                      var.model=c("semi"), var.vars = c(1), knots.v = 3,
                      maxit=1000)
m18[-21]

# semi mean, 3 variance (linear)
m19<-semiVarReg.multi(cc$body_mass_g, x,
                      mean.model="semi", mean.vars = 2, knots.m=3,
                      var.model=c("linear", "linear", "linear"), var.vars = c(1,2,3),knots.v = NULL,
                      maxit=10)
m19[-21]

# semi mean, 3 variance (semi 2 linear 1)
m20<-semiVarReg.multi(cc$body_mass_g, x,
                      mean.model="semi", mean.vars = 2, knots.m=3,
                      var.model=c("semi", "semi", "linear"), var.vars = c(1,2,3),knots.v = c(2,3),
                      maxit=10)
m20[-21]

###########################################################################
# 3 linear mean, constant variance
m21<-semiVarReg.multi(cc$body_mass_g, x,
                      mean.model=c("linear", "linear", "linear"), mean.vars = c(1,2,3),
                      var.model=c("constant"), #var.vars = NULL,
                      maxit=10)
m21[-21]
# 3 linear mean, one variance (linear)
m22<-semiVarReg.multi(cc$body_mass_g, x,
                      mean.model=c("linear", "linear", "linear"), mean.vars = c(1,2,3),
                      var.model=c("linear"), var.vars = c(1),
                      maxit=10)
m22[-21]

# 3 linear mean, one variance (semi)
m23<-semiVarReg.multi(cc$body_mass_g, x,
                      mean.model=c("linear", "linear", "linear"), mean.vars = c(1,2,3),
                      var.model=c("linear"), var.vars = c(1),
                      maxit=10)
m23[-21]

# 3 linear mean, 3 variance (linear)
m24<-semiVarReg.multi(cc$body_mass_g, x,
                      mean.model=c("linear", "linear", "linear"), mean.vars = c(1,2,3),
                      var.model=c("linear", "linear", "linear"), var.vars = c(1,2,3),knots.v = NULL,
                      maxit=10)
m24[-21]

# 3 linear mean, 3 variance (semi 2 linear 1)
m25<-semiVarReg.multi(cc$body_mass_g, x,
                      mean.model=c("linear", "linear", "linear"), mean.vars = c(1,2,3),
                      var.model=c("semi", "semi", "linear"), var.vars = c(1,2,3),knots.v = c(2,3),
                      maxit=10)
m25[-21]

#######################################################################3###
# 1 linear 3 semi mean, constant variance
m26<-semiVarReg.multi(cc$body_mass_g, x,
                      mean.model=c("semi", "linear", "semi"), mean.vars = c(1,2,3),knots.m=c(3,2),
                      var.model=c("constant"), #var.vars = NULL,
                      maxit=10)
m26[-21]

# 1 linear 3 semi mean, one variance (linear)
m27<-semiVarReg.multi(cc$body_mass_g, x,
                      mean.model=c("linear", "semi", "semi"), mean.vars = c(3,2,1),knots.m=c(2,4),
                      var.model=c("linear"), var.vars = c(1),
                      maxit=10)
m27[-21]

# 1 linear 3 semi mean, one variance (semi)
m28<-semiVarReg.multi(cc$body_mass_g, x,
                      mean.model=c("semi", "linear", "semi"), mean.vars = c(4,5,3),knots.m = c(3,4),
                      var.model=c("linear"), var.vars = c(1),
                      maxit=10)
m28[-21]

# 1 linear 3 semi mean, 3 variance (linear)
m29<-semiVarReg.multi(cc$body_mass_g, x,
                      mean.model=c("semi", "linear", "semi"), mean.vars = c(4,5,3),knots.m = c(3,4),
                      var.model=c("linear", "linear", "linear"), var.vars = c(1,2,3),knots.v = NULL,
                      maxit=10)
m29[-21]

# 1 linear 3 semi mean, 3 variance (semi 2 linear 1)
m30<-semiVarReg.multi(cc$body_mass_g, x,
                      mean.model=c("semi", "linear", "semi"), mean.vars = c(4,5,3),knots.m = c(3,4),
                      var.model=c("linear", "semi", "linear"), var.vars = c(1,2,3),knots.v = c(3),
                      maxit=10)
m30[-21]

