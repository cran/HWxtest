## ----setup, include=FALSE------------------------------------------------
library(knitr)
library(knitcitations)
cleanbib()
cite_options(max.names = 1,longnamesfirst=F,  bib.style="authoryear", style="markdown", hyperlink = "to.doc", no.print.fields = c("ISSN", "note"))
library(HWxtest)
library(adegenet)
figw <- 7
figh <- 5
set.seed(60823316) 

## ----refs, include=FALSE-------------------------------------------------
bib <- read.bibtex("bibHW.txt")
engels2009 <- "10.1534/genetics.109.108977"
levene1949 <- "10.1214/aoms/1177730093"
haldane1954 <- "10.1007/BF02985085"
louis1987 <- "10.2307/2531534"
guo1992 <- "10.2307/2532296"
ward2014 <- "10.1093/biostatistics/kxt028"
rousset1995 <- bib[["rousset1995"]]
robertson1984 <- bib[["robertson1984"]]
genepop007 <- bib[["genepop007"]]
adegenet <- "10.1093/bioinformatics/btn129"
pegas <- "10.1093/bioinformatics/btp696"
morin2012 <- bib[["morin2012"]]
citep(morin2012)
gail1977 <- bib[["gail1977"]]
ez1989 <- bib[["ez1989"]]
olsen2014 <- bib[["olsen2014"]]
hart2012 <- bib[["hart2012"]]

## ----set-options, echo=FALSE, cache=FALSE-----------------------------------------------------------------------------
options(width = 120)
library(parallel)
coreCount <- detectCores()
if(coreCount > 1) options(mc.cores=2)
#CRAN seems to reject cores more than 2
#options(mc.cores=detectCores())

## ----include=FALSE----------------------------------------------------------------------------------------------------
obs <- c(83, 49, 18, 74, 34, 21)

## ----entera-----------------------------------------------------------------------------------------------------------
obs <- c(83, 49, 18, 74, 34, 21)

## ----call1------------------------------------------------------------------------------------------------------------
result <- hwx.test(obs)
result

## ----plot1, fig.width=figw, fig.height=figh---------------------------------------------------------------------------
hwx.test(obs, histobins=T)

## ----plot2, fig.width=figw, fig.height=figh---------------------------------------------------------------------------
hwx.test(obs, histobins=T, detail=0, statName="U")

## ----testwhales-------------------------------------------------------------------------------------------------------
data(whales.genind)
wtest <- hwx.test(whales.genind)

## ----dfwhales---------------------------------------------------------------------------------------------------------
dfwhales <- hwdf(wtest)
dfwhales[1:10,]

## ----getcounts--------------------------------------------------------------------------------------------------------
counts1 <- wtest$P1$Bmys42aK46_R225_K232$genotypes
counts1

## ----biggerB----------------------------------------------------------------------------------------------------------
hwx.test(counts1, detail=1, B=1000000)

## ----bigger cutoff----------------------------------------------------------------------------------------------------
counts2 <- wtest$P1$Bmys43Y237_Y377$genotypes
hwx.test(counts2, detail=1, cutoff=2e8)

## ----Udataframe-------------------------------------------------------------------------------------------------------
hwdf(wtest, statName="U")[1:10,]

## ----urchins----------------------------------------------------------------------------------------------------------
urchin.url <- "http://tinyurl.com/ku4fq7m"
urchin.data <- genepop.to.genind(urchin.url)
hwdf(hwx.test(urchin.data))

## ----LLR-vs-asymp, fig.width=figh, fig.height=figh--------------------------------------------------------------------
df <- hwdf(wtest, showAsymptoticX2=T)
pLLR <- df[[1]]
pAsy <- df[[10]]
plot(pLLR, pAsy, xlab="True P value", ylab = "Asymototic approximation")
lines(c(0,1), c(0,1))

## ----zoom, fig.width=figh, fig.height=figh----------------------------------------------------------------------------
lim <- c(0,.2)
plot(pLLR, pAsy, xlim=lim, ylim=lim, xlab="True P value", ylab = "Asymototic approximation")
lines(c(0,1), c(0,1))

## ----useX2, fig.width=figh, fig.height=figh---------------------------------------------------------------------------
dfx <- hwdf(wtest, statName="Chisq")
pX2 <- dfx[[1]]
plot(pX2, pAsy, xlab="True P value", ylab = "Asymototic approximation")
lines(c(0,1), c(0,1))
plot(pX2, pAsy, xlim=lim, ylim=lim, xlab="True P value", ylab = "Asymototic approximation")
lines(c(0,1), c(0,1))

## ----Bmy26------------------------------------------------------------------------------------------------------------
wtest$P2$Bmy26

## ----Bmy26counts, fig.width=figw, fig.height=figh---------------------------------------------------------------------
counts <- wtest$P2$Bmy26$genotypes
hwx.test(counts, detail=0, statName="Chisq", histobins=T, histobounds=c(50, 250), B=1e6)

## ----results="asis", echo=FALSE---------------------------------------------------------------------------------------
bibliography(style=NULL)

