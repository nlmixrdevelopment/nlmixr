par.1cmt.CL <- function(lCL, lV){

CL = exp(lCL)
V = exp(lV)

}

par.1cmt.CL.oral <- function(lCL, lV, lKA){

CL = exp(lCL)
V = exp(lV)
KA = exp(lKA)

}

par.1cmt.CL.oral.tlag <- function(lCL, lV, lKA, lTLAG){

CL = exp(lCL)
V = exp(lV)
KA = exp(lKA)
TLAG = exp(lTLAG)

}

par.2cmt.CL <- function(lCL, lV, lCLD, lVT){

CL = exp(lCL)
V = exp(lV)
CLD = exp(lCLD)
VT = exp(lVT)

}

par.2cmt.CL.oral <- function(lCL, lV, lCLD, lVT, lKA){

CL = exp(lCL)
V = exp(lV)
CLD = exp(lCLD)
VT = exp(lVT)
KA = exp(lKA)

}

par.2cmt.CL.oral.tlag <- function(lCL, lV, lCLD, lVT, lKA, lTLAG){

CL = exp(lCL)
V = exp(lV)
CLD = exp(lCLD)
VT = exp(lVT)
KA = exp(lKA)
TLAG = exp(lTLAG)

}

par.3cmt.CL <- function(lCL, lV, lCLD, lVT, lCLD2, lVT2){

CL = exp(lCL)
V = exp(lV)
CLD = exp(lCLD)
VT = exp(lVT)
CLD2 = exp(lCLD2)
VT2 = exp(lVT2)

}

par.3cmt.CL.oral <- function(lCL, lV, lCLD, lVT, lCLD2, lVT2, lKA){

CL = exp(lCL)
V = exp(lV)
CLD = exp(lCLD)
VT = exp(lVT)
CLD2 = exp(lCLD2)
VT2 = exp(lVT2)
KA = exp(lKA)

}

par.3cmt.CL.oral.tlag <- function(lCL, lV, lCLD, lVT, lCLD2, lVT2, lKA, lTLAG){

CL = exp(lCL)
V = exp(lV)
CLD = exp(lCLD)
VT = exp(lVT)
CLD2 = exp(lCLD2)
VT2 = exp(lVT2)
KA = exp(lKA)
TLAG = exp(lTLAG)

}

par.1cmt.micro <- function(lKE, lV){

KE = exp(lKE)
V = exp(lV)

}

par.1cmt.micro.oral <- function(lKE, lV, lKA){

KE = exp(lKE)
V = exp(lV)
KA = exp(lKA)

}

par.1cmt.micro.oral.tlag <- function(lKE, lV, lKA, lTLAG){

KE = exp(lKE)
V = exp(lV)
KA = exp(lKA)
TLAG = exp(lTLAG)

}

par.2cmt.micro <- function(lKE, lV, lK12, lK21){

KE = exp(lKE)
V = exp(lV)
K12 = exp(lK12)
K21 = exp(lK21)

}

par.2cmt.micro.oral <- function(lKE, lV, lK12, lK21, lKA){

KE = exp(lKE)
V = exp(lV)
K12 = exp(lK12)
K21 = exp(lK21)
KA = exp(lKA)

}

par.2cmt.micro.oral.tlag <- function(lKE, lV, lK12, lK21, lKA, lTLAG){

KE = exp(lKE)
V = exp(lV)
K12 = exp(lK12)
K21 = exp(lK21)
KA = exp(lKA)
TLAG = exp(lTLAG)

}

par.3cmt.micro <- function(lKE, lV, lK12, lK21, lK13, lK31){

KE = exp(lKE)
V = exp(lV)
K12 = exp(lK12)
K21 = exp(lK21)
K13 = exp(lK13)
K31 = exp(lK31)

}

par.3cmt.micro.oral <- function(lKE, lV, lK12, lK21, lK13, lK31, lKA){

KE = exp(lKE)
V = exp(lV)
K12 = exp(lK12)
K21 = exp(lK21)
K13 = exp(lK13)
K31 = exp(lK31)
KA = exp(lKA)

}

par.3cmt.micro.oral.tlag <- function(lKE, lV, lK12, lK21, lK13, lK31, lKA, lTLAG){

KE = exp(lKE)
V = exp(lV)
K12 = exp(lK12)
K21 = exp(lK21)
K13 = exp(lK13)
K31 = exp(lK31)
KA = exp(lKA)
TLAG = exp(lTLAG)

}

.getParfn <- function(oral, ncmt, parameterization, tlag)
{
	x<-sprintf("par.%dcmt.%s%s%s", ncmt,
		c("CL", "micro")[parameterization],
		c("",".oral")[oral+1],
		c("",".tlag")[tlag+1]
	)

	eval(parse(text=x))
}

#get.parfn(oral, ncmt, parameterization, tlag)
