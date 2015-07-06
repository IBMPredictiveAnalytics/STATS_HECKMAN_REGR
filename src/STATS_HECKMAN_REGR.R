#Licensed Materials - Property of IBM
#IBM SPSS Products: Statistics General
#(c) Copyright IBM Corp. 2012, 2014
#US Government Users Restricted Rights - Use, duplication or disclosure 
#restricted by GSA ADP Schedule Contract with IBM Corp.

# Author: JKP, IBM SPSS
# Version = 1.0.1

# history
# 09-Aug-2012 - original version


helptext="The STATS HECKMAN REGR command requires the R Integration Plug-in
and the R sampleSelection package.

This procedure computes Heckman-style censored regression or switching regressions.

STATS HECKMAN REGR
/SELECTION SDEPENDENT=selectionvarname SINDEPENDENT=predictorlist
/OUTCOME ODEPENDENT=outcomevarname OINDEPENDENT=predictorlist
		[OINDEPENDENT2=predictorlist]
[/OPTIONS [METHOD={ML* | TWOSTEP}] [STARTVALUES= list of starting values]
  [DEBUGOUTPUT=number]
	[/SAVE DATASET=datasetname ID=idvariable [SELFIT={YES|NO*} OUTFIT={YES|NO*} SELRES={YES|NO*} OUTRES={YES|NO*}]]
[/HELP]

Example:
STATS HECKMAN REGR SDEPENDENT=employed 
SINDEPENDENT = age education gender
/OUTCOME ODEPENDENT=salary OINDEPENDENT=age education gender experience.

The first form of the command applies when some observations are censored and, hence,
the dependent variable is not observed.  The second form applies when there are two
regimes, and which one applies depends on a set of specified factors.

The SELECTION subcommand specifies the selection model consisting of
a two-valued dependent variable specified as SDEPENDENT and one or more predictors
specified as SINDEPENDENT.  This equation is estimated as a probit model.  The
smaller value of the dependent variable indicates the cases that are censored.
This variable should be numeric.

The OUTCOME subcommand specifies the model that applies for the uncensored
observations, or two models that apply for switching regression.
ODEPENDENT is the dependent variable, and OINDEPENDENT is the
list of independent variables.  Do not code the censored values of the dependent
variable as missing.

The independent variables are assumed to be scale measurement level except
that strings are interpreted as factors.

If ODEPENDENT2 is specified, two regressions are estimated along
with the probit equation.  Both regressions have the same dependent
variable but may have different lists of independent variables.

The models can be estimated by maximum likelihood (ML) or two-step Heckman
regression (TWOSTEP).  Maximum likelihood is the default.  When TWOSTEP
is used, the inverse Mills ratio from the probit equation is included in the
second step in order to address the selection bias.

STARTVALUES can specify starting values for maximum likelihood estimation.
One value must be supplied for each term in the equations plus the rho and sigma
values in the order as shown in the output.

DEBUGOUT specifies the level of debugging output to display.  The default is 0,
which means no debug output.  Larger values produce more information.

The SAVE subcommand specifies what casewise results to save if any.
If used, the DATASET parameter must specify a dataset name that is not
already in use where these results will be written.

SELFIT OUTFIT SELRES and OUTRES can specify saving the selection
probabilities and predicted outcome values and the corresponding residuals.
NOTE: Missing value deletion may cause the number of cases in the results
dataset to be less than the number of input cases.  Cases with missing values
are discarded unless the missing values occur in the outcome equation when
the outcome variable is unobserved.

Since missing values may cause the number of saved cases to differ 
from the input, specify an ID variable to facilitate matching up the saved values 
with the input dataset.

STATS HECKMAN REGR /HELP prints this information and does nothing else.

Details of the statistical models can be found in
Toomet O, Henningsen A (2008b). 'Sample Selection Models in R: Package sampleSelection.'
Journal of Statistical Software, 27(7). URL http://www.jstatsoft.org/v27/i07/
or the equivalent R vignette.
"

heckman <- function(sdep, sindep, odep, oindep, oindep2=NULL,
    method="ml", startvalues=NULL, debugoutput=0, 
		dataset=NULL, selfit=FALSE, outfit=FALSE, selres=FALSE, outres=FALSE, id=NULL) {

    domain<-"STATS_HECKMAN_REGR"
    setuplocalization(domain)
	gtxt <- function(...) {
		return(gettext(...,domain=domain))
	}

	gtxtf <- function(...) {
		return(gettextf(...,domain=domain))
	}
	
    quiet = tryCatch(suppressMessages(library(sampleSelection)), error=function(e) {
        stop(gtxt("The R sampleSelection package is required but could not be loaded."), call. = FALSE)})
    
    if (sdep %in% c(sindep, odep, oindep, oindep2))
				stop(gtxt("The selection  variable may not appear elsewhere as a variable"))
    if (odep %in% c(sdep, sindep, oindep, oindep2))
				stop(gtxt("The outcome variable may not appear elsewhere as a variable"))
		if (is.null(dataset) & any(selfit, outfit, selres, outres))
				stop(gtxt("A dataset name is required if any items are to be saved"))
		if (!is.null(dataset) & !any(selfit, outfit, selres, outres))
				stop(gtxt("A dataset name was specified but no items to save were specified"))
	
	swi = !is.null(oindep2)   # True for switching regression
    
    allvars = union(c(sdep, sindep, odep, oindep), c(oindep2, id))   # expecting order to be preserved
    # categorical variables are not converted because the R package does not handle singularity
    dta<-spssdata.GetDataFromSPSS(allvars, missingValueToNA = TRUE, factorMode = "none")
		if (swi) {
				dta = dta[complete.cases(dta),]  # all cases must be complete for switching regression
		} else {
				# if saving casewise results and there is an id variable, prune missing cases
				# First find all rows with any missings
				# Then restore ok status if the missing occurs in the outcome equation if
				# the case is not selected.
				# Finally re-invalidate the case if the missing occurred in the selection equation.
				if (!is.null(dataset) & !is.null(id)) {
						selmin = min(dta[,1], na.rm=TRUE)  # this defines the unobserved outcome cases
						odeploc = 1 + length(sindep) + 1
						okout = apply(dta[,match(c(odep, oindep), allvars)], 1, function(v) !any(is.na(v)))  # any missing values in outcome eq
						oksel = apply(dta[, c(1:length(sindep)+1)], 1, function(v) !any(is.na(v))) # any missing in sel eq
						okout = ifelse(dta[, 1] == selmin, TRUE, okout)
						dta = dta[okout & oksel,]
				}
		}

    if (method == "twostep") 
				method = "2step"
    rhs = paste(sindep, collapse="+")
    fs = paste(sdep, rhs, sep="~")
    rhs = paste(oindep, collapse="+")
    fo = paste(odep, rhs, sep="~")
		if (swi) {
				rhs = paste(oindep2, collapse="+")
				fo2 = paste(odep, rhs, sep="~")
				modeltype = gtxt("Switching Regression")
		} else {
				modeltype = "Sample Selection"
		}
		if (swi) {
				# selection cannot handle the two formulas as a list variable, so we must
				# write these out explicitly :-(
				res = tryCatch(selection(as.formula(fs), list(as.formula(fo), as.formula(fo2)), data = dta, method=method,
				start=startvalues, print.level=debugoutput), 
        error = function(e) stop(e$message, call. = FALSE))
		} else {
    res = tryCatch(selection(as.formula(fs), as.formula(fo), data = dta, method=method,
				start=startvalues, print.level=debugoutput), 
        error = function(e) stop(e$message, call. = FALSE))
		}
    ressum = summary(res)
    attributes(ressum$estimate)$dimnames[[2]][[1]] = gtxt("Estimate")	
    attributes(ressum$estimate)$dimnames[[2]][[2]] = gtxt("Std. Error")
    attributes(ressum$estimate)$dimnames[[2]][[3]] = gtxt("t Value")
    attributes(ressum$estimate)$dimnames[[2]][[4]] = gtxt("Sig.")
    ress = ressum$estimate[ressum$param$index$betaS,]

		if (swi) {
				if (method == "ml") {
						reso1 = ressum$estimate[ressum$param$index$betaO1,]
						reso2 = ressum$estimate[ressum$param$index$betaO2,]
				} else {
						reso1 = ressum$estimate[c(ressum$param$index$betaO1, ressum$param$index$Mills1),]
						reso2 = ressum$estimate[c(ressum$param$index$betaO2, ressum$param$index$Mills2),]
				}
				R21 = rounder(ressum$rSquared$R21)
				R22 = rounder(ressum$rSquared$R22)
				adjR21 = rounder(ressum$rSquared$R2adj1)
				adjR22 = rounder(ressum$rSquared$R2adj2)
				N1 = ressum$param$N1
				N2 = ressum$param$N2
		} else {
				reso =	ressum$estimate[ressum$param$index$outcome,]
				R2 = rounder(ressum$rSquared$R2)
				adjR2 = rounder(ressum$rSquared$R2adj)
				nCensored = ressum$param$N0
				nObserved = ressum$param$N1
		}
    n = ressum$param$nObs
    nParam = sum(ressum$activePar)
    df = ressum$param$df
		msg = ressum$returnMessage
		if (is.null(msg))
				msg = "---"
    if (method == "ml") {
				likelihood = sprintf("%.4f", ressum$loglik)
				method = "Maximum Likelihood"
				if (swi) {
						sigma1 = ressum$estimate[ressum$param$index$errTerms[1]]
						rho1 = ressum$estimate[ressum$param$index$errTerms[2]]
						sigma2 = ressum$estimate[ressum$param$index$errTerms[3]]
						rho2 = ressum$estimate[ressum$param$index$errTerms[4]]
				}	 else {
						rho = ressum$estimate[ressum$param$index$errTerms[['rho']]]
						if (length(ressum$param$index$errTerms) == 2) {
								sigma = ressum$estimate[ressum$param$index$errTerms[['sigma']]]
						}	else {
						sigma = NA
						}
				}
    } else {  # Heckman estimation
				method = "Two-Step Heckman"
				likelihood = "."
				if (swi) {
					sigma1 = ressum$estimate[ressum$param$index$sigma1]
					sigma2 = ressum$estimate[ressum$param$index$sigma2]
					rho1 = ressum$estimate[ressum$param$index$rho1]
					rho2 = ressum$estimate[ressum$param$index$rho2]	
				} else {
					sigma = ressum$estimate[ressum$param$index$sigma]
					rho = ressum$estimate[ressum$param$index$rho]
				}
    }
		
		# result saving
		if (!is.null(dataset)) {   # saving something
				item = 0
				dictlist = list()
				ds = data.frame(c(1:ressum$param$nObs))

				if (!is.null(id)) {  # copy id variable with same meta information to new dataset
						item = item + 1
						iddic = spssdictionary.GetDictionaryFromSPSS(id)
						dictlist[[item]] = c("ID", iddic[2,1], as.integer(iddic[3,1]), iddic[4,1], iddic[5,1])
						ds[item] = dta[,length(allvars)]  # id variable is last
				}
				if (selfit) {
						item = item + 1
						dictlist[[item]] = c("selection_fitted", gtxt("Selection Fitted Values"), 0, "F8.2", "scale")
						ds[item] = fitted(res, part="selection")
				}
				if (outfit) {
						item = item + 1
						dictlist[[item]] = c("outcome_fitted", gtxt("Outcome Fitted Values"), 0, "F8.2", "scale")
						ds[item] = fitted(res, part="outcome")
				}
				if (selres) {
						item = item + 1
						dictlist[[item]] = c("selection_residuals", gtxt("Selection Residuals"), 0, "F8.2", "scale")
						ds[item] = residuals(res, part="selection")
				}
				if (outres) {
						item = item + 1
						dictlist[[item]] = c("outcome_residuals", gtxt("Outcome Residuals"), 0, "F8.2", "scale")
						ds[item] = residuals(res, part="outcome")
				}

				dict = spssdictionary.CreateSPSSDictionary(dictlist)
				tryCatch({spssdictionary.SetDictionaryToSPSS(dataset, dict)
						spssdata.SetDataToSPSS(dataset, ds)}, 
						finally=spssdictionary.EndDataStep())
		}

        
    # print results
    # 
    StartProcedure(gtxt("Heckman Regression"), "STATS HECKMAN REGR") 

    # summary statistics
		scaption = gtxt("Implemented using package sampleSelection by O. Toomet and A.Henningsen")
		if (swi) {
				lbls = c(gtxt("Selection Variable"), gtxt("Outcome Variable"), gtxt("Model Type"), gtxt("Number of Valid Cases"), 
						gtxt("Number of Outcome 1 Cases"), gtxt("Number of Outcome 2 Cases"), gtxt("Number of Parameters"),
						gtxt("D.F"), gtxt("R-Squared 1"), gtxt("Adjusted R-Squared 1"), 
						gtxt("R-Squared 2"), gtxt("Adjusted R-Squared 2"),
						gtxt("Method"), gtxt("Computations"), gtxt("Dataset Created"))
				vals = c(sdep, odep, modeltype, n, N1, N2, nParam, df, R21, adjR21, R22, adjR22, method, 
				msg, ifelse(is.null(dataset), "---", dataset))
				spsspivottable.Display(data.frame(vals, row.names=lbls), title = gtxt("Summary Statistics"),
				collabels=c(gtxt("Summary")), templateName="HECKMANFIT", outline=gtxt("Summary"),
				caption = scaption)
				
				caption = sprintf(gtxt("Selection Variable: %s"), sdep)
				spsspivottable.Display(ress, title=gtxt("Probit Selection Estimates"), templateName="HECKMANSEL",
						outline=gtxt("Selection Estimates"), caption=caption)
			
				caption = sprintf(gtxt("Outcome Variable: %s, Sigma1: %.4f, Rho1: %.4f"), odep, sigma1, rho1)    
				spsspivottable.Display(reso1, title=gtxt("Outcome 1 Estimates"), templateName="HECKMANOUT",
			outline=gtxt("Outcome 1 Estimates"), caption=caption)
			
				caption = sprintf(gtxt("Outcome Variable: %s, Sigma2: %.4f, Rho2: %.4f"), odep, sigma2, rho2)    
				spsspivottable.Display(reso2, title=gtxt("Outcome 2 Estimates"), templateName="HECKMANOUT2",
			outline=gtxt("Outcome 2 Estimates"), caption=caption)
		} else {
		    lbls = c(gtxt("Selection Variable"), gtxt("Outcome Variable"), gtxt("Model Type"), gtxt("Number of Valid Cases"), 
						gtxt("Number of Censored Cases"), gtxt("Number of Observed Cases"), gtxt("Number of Parameters"),
						gtxt("D.F"), gtxt("R-Squared"), gtxt("Adjusted R-Squared"), gtxt("Method"), gtxt("Log Likelihood"),
						gtxt("Computations"), gtxt("Dataset Created"))
        
				vals = c(sdep, odep, modeltype, n, nCensored, nObserved, nParam, df, R2, adjR2, method, 
						likelihood, msg, ifelse(is.null(dataset), "---", dataset))
				spsspivottable.Display(data.frame(vals, row.names=lbls), title = gtxt("Summary Statistics"),
						collabels=c(gtxt("Summary")), templateName="HECKMANFIT", outline=gtxt("Summary"),
						caption = scaption)
				
				caption = sprintf(gtxt("Selection Variable: %s"), sdep)
				spsspivottable.Display(ress, title=gtxt("Probit Selection Estimates"), templateName="HECKMANSEL",
						outline=gtxt("Selection Estimates"), caption=caption)
			
				caption = sprintf(gtxt("Outcome Variable: %s, Sigma: %.4f, Rho: %.4f"), odep, sigma, rho)    
				spsspivottable.Display(reso, title=gtxt("Outcome Estimates"), templateName="HECKMANOUT",
						outline=gtxt("Outcome Estimates"), caption=caption)
	}
    spsspkg.EndProcedure()

    # clean up workspace
    res <- tryCatch(rm(list=ls()), warning = function(e) {return(NULL)})
}

rounder <- function(s) {
		if (!is.null(s)) {
				s = round(s, 4)
		} else {
				s = "."
		}
		return(s)
}
		
# override for api to account for extra parameter in V19 and beyond
StartProcedure <- function(procname, omsid) {
    if (substr(spsspkg.GetSPSSVersion(),1, 2) >= 19) {
       spsspkg.StartProcedure(procname, omsid)
    }
    else {
       spsspkg.StartProcedure(omsid)
    }
}
setuplocalization = function(domain) {
    # find and bind translation file names
    # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
    # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

    fpath = Find(file.exists, file.path(.libPaths(), paste(domain, ".R", sep="")))
    bindtextdomain(domain, file.path(dirname(fpath), domain, "lang"))
} 

Run<-function(args){
    
    cmdname = args[[1]]
    args <- args[[2]]
    oobj<-spsspkg.Syntax(templ=list(
            spsspkg.Template("SDEPENDENT", subc="SELECTION",  ktype="existingvarlist", var="sdep", islist=FALSE),
            spsspkg.Template("SINDEPENDENT", subc="SELECTION",  ktype="existingvarlist", var="sindep", islist=TRUE),
            spsspkg.Template("ODEPENDENT", subc="OUTCOME",  ktype="existingvarlist", var="odep", islist=FALSE),
            spsspkg.Template("OINDEPENDENT", subc="OUTCOME",  ktype="existingvarlist", var="oindep", 
                islist=TRUE),
						spsspkg.Template("OINDEPENDENT2", subc="OUTCOME",  ktype="existingvarlist", var="oindep2", 
                islist=TRUE),
            spsspkg.Template("METHOD", subc="OPTIONS", ktype="str", var="method", 
                vallist=list("ml", "twostep")),
						spsspkg.Template("DATASET", subc="SAVE", ktype="varname", var="dataset"),
						spsspkg.Template("SELFIT", subc="SAVE", ktype="bool", var="selfit"),		
						spsspkg.Template("OUTFIT", subc="SAVE", ktype="bool", var="outfit"),
						spsspkg.Template("SELRES", subc="SAVE", ktype="bool", var="selres"),
						spsspkg.Template("OUTRES", subc="SAVE", ktype="bool", var="outres"),
						spsspkg.Template("ID", subc="SAVE", ktype="existingvarlist", var="id", islist=FALSE),
            spsspkg.Template("STARTVALUES", subc="OPTIONS", ktype="float", var="startvalues", islist=TRUE),
            spsspkg.Template("DEBUGOUTPUT", subc="OPTIONS", ktype="int", var="debugoutput")
                ))        
    if ("HELP" %in% attr(args,"names"))
        #writeLines(helptext)
        helper(cmdname)
    else
        res <- spsspkg.processcmd(oobj,args,"heckman")
}

helper = function(cmdname) {
    # find the html help file and display in the default browser
    # cmdname may have blanks that need to be converted to _ to match the file
    
    fn = gsub(" ", "_", cmdname, fixed=TRUE)
    thefile = Find(file.exists, file.path(.libPaths(), fn, "markdown.html"))
    if (is.null(thefile)) {
        print("Help file not found")
    } else {
        browseURL(paste("file://", thefile, sep=""))
    }
}
if (exists("spsspkg.helper")) {
assign("helper", spsspkg.helper)
}
