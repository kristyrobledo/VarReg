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
# constant location, constant scale2, constant shape
m1<-lssVarReg.multi(y, x,
                    locationmodel="constant", # location.vars = c(1),
                      scale2model="constant",   #  scale2.vars = c(1),
                    shapemodel="constant")    #  shape.vars = c(1),
                    #knots.l=2, knots.sc=2, knots.sh, degree=2,
m1[-21]

# constant location, constant scale2, linear shape
m2<-lssVarReg.multi(y, x,
                    locationmodel="constant", # location.vars = c(1),
                    scale2model="constant",   #  scale2.vars = c(1),
                    shapemodel="linear", shape.vars = 3)    #  shape.vars = c(1),
#knots.l=2, knots.sc=2, knots.sh, degree=2,
m2[-21]

# constant location, constant scale2, semi shape
m3<-lssVarReg.multi(y, x,
                    locationmodel="constant", # location.vars = c(1),
                    scale2model="constant",   #  scale2.vars = c(1),
                    shapemodel="semi", shape.vars = 3, knots.sh = 2, degree = 2)    #  shape.vars = c(1),
#knots.l=2, knots.sc=2, knots.sh, degree=2,
m3[-21]

# constant location, constant scale2, 2 linear shape
m2<-lssVarReg.multi(y, x,
                    locationmodel="constant", # location.vars = c(1),
                    scale2model="constant",   #  scale2.vars = c(1),
                    shapemodel=c("linear","linear"), shape.vars = c(3,1))    #  shape.vars = c(1),
#knots.l=2, knots.sc=2, knots.sh, degree=2,
m2[-21]
# constant location, constant scale2, linear and semi shape
m2<-lssVarReg.multi(y, x,
                    locationmodel="constant", # location.vars = c(1),
                    scale2model="constant",   #  scale2.vars = c(1),
                    shapemodel=c("linear", "semi"), shape.vars = c(2,3), knots.sh = 1 )    #  shape.vars = c(1),
#knots.l=2, knots.sc=2, knots.sh, degree=2,
m2[-21]

# constant location, constant scale2, 2 semi shape
m2<-lssVarReg.multi(y, x,
                    locationmodel="constant", # location.vars = c(1),
                    scale2model="constant",   #  scale2.vars = c(1),
                    shapemodel=c("semi", "semi"), shape.vars = c(2,3), knots.sh = c(1,2) )

m2[-21]

# constant location, linear scale2, constant shape
m2<-lssVarReg.multi(y, x,
                    locationmodel="constant", # location.vars = c(1),
                    scale2model="linear",  scale2.vars = c(1),
                    shapemodel=c("constant") )

# constant location, linear scale2, linear shape

m2<-lssVarReg.multi(y, x,
                    locationmodel="constant", # location.vars = c(1),
                    scale2model="linear",  scale2.vars = c(1),
                    shapemodel=c("linear"),shape.vars = c(2) )


m2<-lssVarReg.multi(y, x,
                    locationmodel="constant", # location.vars = c(1),
                    scale2model="linear",  scale2.vars = c(1),
                    shapemodel=c("semi", "semi"), shape.vars = c(2,3), knots.sh = c(1,2) )


m2<-lssVarReg.multi(y, x,
                    locationmodel="constant", # location.vars = c(1),
                    scale2model="semi",  scale2.vars = c(1), knots.sc=2,
                    shapemodel=c("semi", "semi"), shape.vars = c(2,3), knots.sh = c(1,2) )


m2[-21]


m2<-lssVarReg.multi(y, x,
                    locationmodel="constant", # location.vars = c(1),
                    scale2model=c("linear", "semi"),  scale2.vars = c(1,2), knots.sc=2,
                    shapemodel=c("linear"), shape.vars = c(2), print.it = TRUE)

m2[-21]
## runs but estiamtes are crappy


##this is failing....

m2<-lssVarReg.multi(y, x,
                    locationmodel=c("linear", "semi"),  location.vars = c(1,2), knots.l=2, # location.vars = c(1),
                    scale2model=c("linear"), scale2.vars = c(1),
                    shapemodel=c("linear"), shape.vars = c(2), print.it = TRUE,
                    scale2.init = c(2,2),
                    shape.init = c(3,3))

m2[-21]
## need to build in the initial estimates option! ->DONE



m2<-lssVarReg.multi(y, x,
                    locationmodel="constant",   # location.vars = c(1),
                    scale2model="constant",
                    shapemodel=c("linear", "semi"), shape.vars = c(2,3), knots.sh = 2, print.it = TRUE)

m2[-21]

# constant location, linear scale2, semi shape
# constant location, linear scale2, 2 linear shape
# constant location, linear scale2, linear and semi shape
# constant location, linear scale2, 2 semi shape

# constant location, semi scale2, constant shape
# constant location, semi scale2, linear shape
# constant location, semi scale2, semi shape
# constant location, semi scale2, 2 linear shape
# constant location, semi scale2, linear and semi shape
# constant location, semi scale2, 2 semi shape

# constant location, semi scale2, constant shape
# constant location, semi scale2, linear shape
# constant location, semi scale2, semi shape
# constant location, semi scale2, 2 linear shape
# constant location, semi scale2, linear and semi shape
# constant location, semi scale2, 2 semi shape



allshapemodels <-list("constant",
                      "linear",
                      "semi",
                      c("linear", "linear"),
                      c("semi", "linear"),
                      c("semi", "semi"))
allshapevars<-list(NULL,
                   "2",
                   "2",
                   c(3,2),
                   c(3,4),
                   c(4,2))

allshapeknots<-list(NULL,
                    NULL,
                    3,
                    NULL,
                    2,
                    c(1,2))

fmodels<-list()


## constant location, linear scale2, varying shapes:

for (i in 1:6){
  print(i)
  fmodels[[i]]<-lssVarReg.multi(y, x,
                    locationmodel="constant", # location.vars = c(1),
                    scale2model="constant",   #  scale2.vars = c(1),
                    shapemodel=allshapemodels[i],
                    shape.vars =allshapevars[i],
                    knots.sh =allshapeknots[i] )
}

fmodels[[2]]<-lssVarReg.multi(y, x,
                            locationmodel="constant", # location.vars = c(1),
                            scale2model="constant",   #  scale2.vars = c(1),
                            shapemodel=allshapemodels[2],
                            shape.vars =allshapevars[2],
                            knots.sh =allshapeknots[2] )


