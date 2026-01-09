################################################################################
### Function that calculates the probability of 2-year mortality for patients ##
### treated with dialysis or conservative care #################################
################################################################################
prob_D_and_CC <- function(X, arm="D"){
  if (arm=="D"){
    # extended dialysis model
    a <- -4.6099
    b <- c(0.0737,  # age
           0.0214,  # eGFR
           -0.1185, # malignancy present
           0.0470,  # diabetes mellitus
           0.3588,  # ischaemic heart disease
           1.1720,  # left ventricular dysfunction
           1.0187,  # peripheral vascular disease
           -0.0245, # gender
           -0.0820, # serum albumin
           0.0015)  # CRP
  } else if (arm=="CC"){
    # extended conservative care model
    a <- 1.5472
    b <- c(0.0374, # age
            -0.1049, # eGFR
            0.4742,  # malignancy present
            -0.1753, # diabetes mellitus
            -0.0843, # ischaemic heart disease
            1.2749,  # left ventricular dysfunction
            0.7225,  # peripheral vascular disease
            -0.1839, # gender
            -0.0927, # serum albumin
            0.0012)  # CRP
  }
  
  # prognostic index
  PI <- a + b%*%t(X[,-1])
  
  return(plogis(as.numeric(PI)))
}
# patient characteristics
X <- data.frame(A=c(1, 0),         # treatment assignment
                age=c(70, 70),     # age
                eGFR=c(10, 10),    # eGFR
                malignancy=c(1, 1),# malignancy present
                diabetes=c(0, 0),  # diabetes mellitus
                ihd=c(0, 0),       # ischaemic heart disease
                lvd=c(0, 0),       # left ventricular dysfunction
                pvd=c(0, 0),       # peripheral vascular disease
                sex=c(1, 1),       # gender (female)
                albu=c(40, 40),     # serum albumin
                CRP=c(50, 50))     # CRP

# obtain probabilities for identical patient with dialysis or conservative care
prob_D_and_CC(X, arm="D")*100