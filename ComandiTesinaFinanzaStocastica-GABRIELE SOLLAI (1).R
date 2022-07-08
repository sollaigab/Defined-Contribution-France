library(StMoMo)
library(demography)
library(lifecontingencies)

#Caricamento dati
fr=read.demogdata("franciamx.txt", "franciaex.txt", type="mortality", label="France")
fr
#Plot delle serie
par(mfrow=c(1,3))
plot(fr,series="male",datatype="rate",  plot.type="time", main="Male rates", xlab="years")
plot(fr,series="female",datatype="rate", plot.type="time", main="Female rates", xlab="years")
plot(fr,"total",datatype="rate", plot.type="time", main="Total rates", xlab="years")
par(mfrow=c(1,3))

#Uso di StMoMo
FR=extract.years(fr,1963:2017)
FRA=extract.ages(FR,0:100)
France=StMoMoData(FRA, series="total")
FranceIni=central2initial(France)
ages.fit=0:100
wxt=genWeightMat(ages=ages.fit, years=FranceIni$years, clip=3)

#Modello LC
LC=lc("logit")
LCfit=fit(LC, data=FranceIni, ages.fit=ages.fit, wxt=wxt)
plot(LCfit)
LCres=residuals(LCfit)
plot(LCres, type="colourmap")
AIC(LCfit)
BIC(LCfit)
LCfit$deviance

#Modello APC
APC=apc("logit")
APCfit=fit(APC, data=FranceIni, ages.fit=ages.fit, wxt=wxt)
plot(APCfit)
APCres=residuals(APCfit)
plot(APCres, type="colourmap")
AIC(APCfit)
BIC(APCfit)
APCfit$deviance

#Modello RH
RH=rh("logit")
RHfit=fit(RH, data=FranceIni, ages.fit=ages.fit, wxt=wxt)
plot(RHfit)
RHres=residuals(RHfit)
plot(RHres, type="colourmap")
AIC(RHfit)
BIC(RHfit)
RHfit$deviance

#Modello CBD
CBD=cbd("logit")
CBDfit=fit(CBD, data=FranceIni, ages.fit=ages.fit, wxt=wxt)
plot(CBDfit)
CBDres=residuals(CBDfit)
plot(CBDres,"scatter" )
plot(CBDres, type="colourmap")
AIC(CBDfit)
BIC(CBDfit)
CBDfit$deviance

#Modello M7
M7=m7("logit")
M7fit=fit(M7, data=FranceIni, ages.fit=ages.fit, wxt=wxt)
plot(M7fit)
M7res=residuals(M7fit)
plot(M7res)
plot(M7res, type="colourmap")
AIC(M7fit)
BIC(M7fit)
M7fit$deviance

#Modello PLAT
f2 <- function(x, ages) mean(ages) - x
constPlat <- function(ax, bx, kt, b0x, gc, wxt, ages){
nYears <- dim(wxt)[2]
x <- ages
t <- 1:nYears
c <- (1 - tail(ages, 1)):(nYears - ages[1])
xbar <- mean(x)
phiReg <- lm(gc ~ 1 + c + I(c ^ 2), na.action = na.omit)
phi <- coef(phiReg)
gc <- gc - phi[1] - phi[2] * c - phi[3] * c ^ 2
kt[2, ] <- kt[2, ] + 2 * phi[3] * t
kt[1, ] <- kt[1, ] + phi[2] * t + phi[3] * (t ^ 2 - 2 * xbar * t)
ax <- ax + phi[1] - phi[2] * x + phi[3] * x ^ 2
ci <- rowMeans(kt, na.rm = TRUE)
ax <- ax + ci[1] + ci[2] * (xbar - x)
kt[1, ] <- kt[1, ] - ci[1]
kt[2, ] <- kt[2, ] - ci[2]
list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)
}
PLAT <- StMoMo(link = "logit", staticAgeFun = TRUE,
periodAgeFun = c("1", f2), cohortAgeFun = "1", constFun = constPlat)
PLATfit <- fit(PLAT, data = FranceIni, ages.fit = ages.fit, wxt = wxt)
plot(PLATfit)
PLATres=residuals(PLATfit)
plot(PLATres)
plot(PLATres, type="colourmap")
AIC(PLATfit)
BIC(PLATfit)
PLATfit$deviance

#Forecasting con LC
LCfT=forecast(LCfit, h=50)
plot(LCfT)

#Forecasting con PLAT
PLATfT=forecast(PLATfit, h=50)
plot(PLATfT)

#Creazione di simulazioni
nsims=100
LCsim=simulate(LCfit, nsim = nsims, h = 50)
PLATsim=simulate(PLATfit, nsim = nsims, h = 50)

#Plot Simulazioni LC
plot(LCfit$years, LCfit$kt[1, ], xlim = range(LCfit$years, LCsim$kt.s$years),
ylim = range(LCfit$kt, LCsim$kt.s$sim[1, , 1:20]), type = "l",
xlab = "year", ylab = "kt", main = "LC model simulations")
matlines(LCsim$kt.s$years, LCsim$kt.s$sim[1, , 1:20], type = "l", lty=1)

#Plot Simulazioni PLAT
plot(PLATfit$years, PLATfit$kt[1, ], xlim = range(PLATfit$years, PLATsim$kt.s$years),
ylim = range(PLATfit$kt, PLATsim$kt.s$sim[1, , 1:20]), type = "l",
xlab = "year", ylab = "kt", main = "PLAT model simulations")
matlines(PLATsim$kt.s$years, PLATsim$kt.s$sim[1, , 1:20], type = "l", lty=1)

#Grafico per una specifica coorte con StMoMo con LC
chosen_cohort=1963
plot(0:54, extractCohort(fitted(LCfit, type = "rates"), cohort = chosen_cohort),
type = "l", log = "y", xlab = "age", ylab = "q(x)",
main = "LC Cohort 1963 mortality rate", xlim = c(0,100), ylim = c(0.00035, 0.3))
#adding fitted projections
lines(55:100, extractCohort(LCfT$rates, cohort = chosen_cohort), lty = 2, lwd=2, col="red")

#Grafico per una specifica coorte con StMoMo con PLAT
chosen_cohort=1963
plot(0:54, extractCohort(fitted(PLATfit, type = "rates"), cohort = chosen_cohort),
type = "l", log = "y", xlab = "age", ylab = "q(x)",
main = " PLAT Cohort 1963 mortality rate", xlim = c(0,100), ylim = c(0.00035, 0.3))
#adding fitted projections
lines(55:100, extractCohort(PLATfT$rates, cohort = chosen_cohort), lty = 2, lwd=2, col="red")

#Mortality rates LC
LC_historical_rates <- extractCohort(fitted(LCfit, type = "rates"),
cohort = chosen_cohort)
LC_forecasted_rates <- extractCohort(LCfT$rates,
cohort = chosen_cohort)
LC_rates_1963 <- c(LC_historical_rates,LC_forecasted_rates)

#Mortality rates PLAT
PLAT_historical_rates <- extractCohort(fitted(PLATfit, type = "rates"),
cohort = chosen_cohort)
PLAT_forecasted_rates <- extractCohort(PLATfT$rates,
cohort = chosen_cohort)
PLAT_rates_1963 <- c(PLAT_historical_rates,PLAT_forecasted_rates)


#Trasformare tassi di mortalità in probabilità di morte
lc_qx_1963<-mx2qx(lc_rates_1963)
PLAT_qx_1963<-mx2qx(PLAT_rates_1963)

#Lifetable con LC
lc_lifetable_1963<-probs2lifetable(probs=lc_qx_1963,type = "qx",
name = paste("LC","1963","lt",sep="_"))
PLAT_lifetable_1963<-probs2lifetable(probs=PLAT_qx_1963,type = "qx",
name = paste("LC","1963","lt",sep="_"))

#Aspettativa di vita esempio
exn(lc_lifetable_1963,x=57)
exn(PLAT_lifetable_1963,x=57)

#Creazione Actuarial Table
lc_acttbl_1963<-new("actuarialtable",x=lc_lifetable_1963@x,lx=lc_lifetable_1963@lx,
interest=0.014,name="LC Actuarial Table")
PLAT_acttbl_1963<-new("actuarialtable",x=PLAT_lifetable_1963@x,lx=PLAT_lifetable_1963@lx,
interest=0.014,name="PLAT Actuarial Table")

#Applicazione metodo IEAM ai modelli
IEAM <- function(acttableAccPeriod, x, beta, i = actuarialtable@interest, j, 
    t, k = 1, payment = "advance", acttablePaymPeriod, i2, delta = 0, type = 0) {
    
    out <- numeric(1)
    if (missing(acttableAccPeriod)) 
        stop("Error! Need an actuarial actuarialtable")
    if (missing(acttablePaymPeriod)) 
        acttablePaymPeriod = acttableAccPeriod
    if (missing(i2)) 
        i2 = i
    if (missing(x)) 
        stop("Error! Need age!")
    if (missing(beta)) 
        stop("Error! Retirement age!")
    if (x > getOmega(acttableAccPeriod)) {
        out = 0
        return(out)
    }
    if (missing(t)) 
        stop("Error! Need t")
    if (missing(j)) 
        stop("Error! Need average salary increase rate")
    if (any(x < 0, beta < 0, t < 0)) 
        stop("Error! Negative parameters")
    if (type == 0) {
        out = (Exn(acttableAccPeriod, x, beta - x, i = i) * (beta - x)/t * axn(acttablePaymPeriod, 
            beta, i = (1 + i2)/(1 + delta) - 1, k = 1) * (1 + j)^(beta - x - 
            1))/(axn(acttablePaymPeriod, x, beta - x, i = (1 + i)/(1 + j) - 
            1, k = 1, payment = "advance"))
    } else {
        out = ((Exn(acttableAccPeriod, x, beta - x, i = i) * (beta - x)/t * 
            axn(acttablePaymPeriod, beta, i = (1 + i2)/(1 + delta) - 1, k = 1) * 
            (1 + j)^(beta - x - 1))/(axn(acttablePaymPeriod, x, beta - x, i, 
            k = 1, payment = "advance")))/(1 + j)^seq(0, beta - x - 1, 1)
    }
    return(out)
}
beta=62 # Beta Retirement age
x=25 #x Age of the insured.
i=0.08 # Interest Rate
t=37 #1/t is the % of the salary, recognized as retirement pension, for each year of service
j=0.01 #  average salary increases (for both growth in wages and promotional salary for seniority)
delta=0.03 #Increase of retirement pension

#Calcolo contribution rate
CRlc=IEAM(lc_acttbl_1963, x, beta, i, j, t, k, delta = 0.03, type = 0)
CRplat=IEAM(PLAT_acttbl_1963, x, beta, i, j, t, k, delta = 0.03, type = 0)

#Calcolo Attività del fondo
aumento=1.01
x=50000*aumento
y=x*aumento
r=y*aumento
e=r*aumento
w=e*aumento
q=w*aumento
t=q*aumento
u=t*aumento
i=u*aumento
a=i*aumento
o=a*aumento
p=o*aumento
s=p*aumento
d=s*aumento
f=d*aumento
g=f*aumento
h=g*aumento
j=h*aumento
k=j*aumento
l=k*aumento
z=l*aumento
v=z*aumento
b=v*aumento
n=b*aumento
m=n*aumento
Q=m*aumento
W=Q*aumento
E=W*aumento
R=E*aumento
T=R*aumento
Y=T*aumento
U=Y*aumento
I=U*aumento
O=I*aumento
P=O*aumento
A=P*aumento
salario=c(50000,x,y,r,e,w,q,t,u,i,a,o,p,s,d,f,g,h,j,k,l,z,v,b,n,m,Q,W,E,R,T,Y,U,I,O,P,A)
salario
salario_tot=salario*100
attività_lc=salario_tot*CRlc
attività_plat=salario_tot*CRplat
plot(attività_plat,type="l",main="Confronto tra Contributi con LC e PLAT",xlab="Numero contributi",ylab="Contributi",col="red")
lines(attività_lc,col="blue")
legend("topleft" , c("Contributi Plat","Contributi LC"), cex=0.8,col=c("red","blue"),lty=1)
write.table(attività_lc, file="attività.csv", quote=F, sep=";", dec=",", na="", row.names=TRUE, col.names=TRUE)
write.table(attività_plat, file="attivitàplat.csv", quote=F, sep=";", dec=",", na="", row.names=TRUE, col.names=TRUE)
