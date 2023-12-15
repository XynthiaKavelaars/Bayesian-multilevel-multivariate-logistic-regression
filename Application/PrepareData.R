# Select variables
SelectedVariables <- Dataset[,c("AGE", "RSBP", "RDELAY", "RXASP", "RXHEP", "STRK14", "H14", "OCCODE")]
SelectedVariables$Trt <- interaction(SelectedVariables$RXASP,SelectedVariables$RXHEP)

# Select cases
# OCCODE: 1 = DEAD, 2 = DEPENDENT, 3 = NOT DEPENDENT, 4 = RECOVERED
# TRT: Y.N = YES ASPIRIN/NO HEPARIN; N.N = NO ASPIRIN/NO HEPARIN
Data <- SelectedVariables[SelectedVariables$Trt %in% c("Y.M", "Y.H", "Y.N") & 
                            SelectedVariables$OCCODE %in% c(2:4),]

# Recode treatment and outcome variables
Data$TRT <- Data$Trt %in% c("Y.M", "Y.H")       # 1 = Aps + Hep, 0 = Asp
Data$DEP6 <- Data$OCCODE == 2                   # 1 = Dependent; 0 = Not dependent

# Sample size
n <- nrow(Data)                                 # Selected subjects
nE <- sum(Data$TRT == 1)                        # Selected subjects in treatment 1
nC <- n - nE                                    # Selected subjects in treatment 0

# Covariate data 
bp <- Data$RSBP
Data$RSBP_C <- Data$RSBP - mean(Data$RSBP)
#XBp <- matrix(seq(min(Data$RSBP),max(Data$RSBP), length.out=1e2), ncol = 1)
#XBp_C <- matrix(seq(min(Data$RSBP_C),max(Data$RSBP_C), length.out=1e2), ncol = 1)

trts <- Data[,c("TRT")]

# Summary of distributions covariate data 
#xBp <- matrix(min(Data$RSBP):max(Data$RSBP), ncol = 1)
#XEBp <- cbind(rep(1,nrow(xBp)), rep(1,nrow(xBp)), xBp, xBp)
#XCBp <- cbind(rep(1,nrow(xBp)), rep(0,nrow(xBp)), xBp, rep(0, nrow(xBp)))

meanBp <- mean(Data$RSBP)
sdBp <- sd(Data$RSBP)

meanBp_E <- mean(Data[Data$TRT == 1, "RSBP"])
meanBp_C <- mean(Data[Data$TRT == 0, "RSBP"])
sdBp_E <- sd(Data[Data$TRT == 1, "RSBP"])
sdBp_C <- sd(Data[Data$TRT == 0, "RSBP"]) 

# Response data
# Reverse-coded! 0 = NO RECURRENT STROKE/NOT DEPENDENT; 1 = RECURRENT STROKE/DEPENDENT
Y <- Data[,c("STRK14", "DEP6")]
Y12 <- matrix(NA, nrow=nrow(Y), ncol=4)
Y12[,1] <- Y[,"STRK14"] == 1 & Y[,"DEP6"] == 1    # Y11
Y12[,2] <- Y[,"STRK14"] == 1 & Y[,"DEP6"] == 0    # Y10
Y12[,3] <- Y[,"STRK14"] == 0 & Y[,"DEP6"] == 1    # Y01
Y12[,4] <- Y[,"STRK14"] == 0 & Y[,"DEP6"] == 0    # Y00
Y12<- Y12 * 1

colnames(Y12) <- c("Yes strk/yes dep", "Yes strk/ No dep", "No strk/Yes dep", "No strk/no dep")
