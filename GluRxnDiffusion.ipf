//IntegrateODE/m=0/u=500 GluatamateKinetic, kk, glukin
// SetDimLabel 1, 0, A, CalKin // Ca  Sets dimension labels to reactant names, see notebook for key
// CalKin[0][%A] = Carest
menu "macros"
  "D3dGaussT"
  "SolveGluRt"
  "D3dGaussX"
  "SetM"
  "BoundFractionVsDist"
  "SetConstants"
end

Function/wave SetConstants()
if(!exists("wBeta"))
make /n=15 wBeta = {0.117,0.423,173,1,1,1.1e6,190,0,1,0,2,10,0,2^8,5e-3}
SetDimLabel 0, 0, SigmaXY, wBeta; SetDimLabel 0, 1, SigmaZ, wBeta; SetDimLabel 0, 2, D, wBeta; SetDimLabel 0, 3, R, wBeta; SetDimLabel 0, 4, M, wBeta
SetDimLabel 0, 5, Kf, wBeta; SetDimLabel 0, 6, Kr, wBeta; SetDimLabel 0, 7, tt, wBeta;SetDimLabel 0, 8, PeakGlu, wBeta
SetDimLabel 0, 9, Xmin, wBeta; SetDimLabel 0, 10, Xmax, wBeta; SetDimLabel 0, 11, nXsteps, wBeta
SetDimLabel 0, 12, Kdecay, wBeta
SetDimLabel 0, 13, NumPntsReceptorR, wBeta
SetDimLabel 0, 14, TimeScaleReceptorR, wBeta
endif

//allow user interaction to edit parameters in table
dowindow /f Parameters
if(!v_flag)
edit /n=Parameters wBeta.l,wBeta
NewPanel/K=2 /n=PauseForUser0 as "Pause for user"; AutoPositionWindow/M=1/R=Parameters
DrawText 21,20,"Edit Parameters in table";
Button button0,pos={80,58},size={92,20},title="Continue"; Button button0,proc=UserContinue
PauseForUser PauseForUser0, Parameters
endif
return wBeta
end

Function UserContinue(ctrlName) : ButtonControl
	String ctrlName
	DoWindow/K PauseForUser0
End

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//SetM sets the constant M such that the peak concentration at x = 0 is PeakGlu
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
function/wave SetM()
wave wBeta = SetConstants()
wBeta[%tt] = 0
wave Cx = D3dGaussX()
wBeta[%M] = wBeta[%PeakGlu]*wBeta[%M]/wavemax(Cx)
D3dGaussX()
return wBeta
end

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//D3dGaussX returns the analytical solution to diffusion in 3D with r varied at fixed t
//with gaussian initial conditions
// initial gaussian is assumed to be symmetric in xy
// the solution assumes z = 0
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
function/wave D3dGaussX()
wave wBeta = SetConstants()
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//Data Dictionary
// Variable vSigmaXY = wBeta[%SigmaXY] //gaussian width in xy in µm
// Variable vSigmaZ = wBeta[%SigmaZ] //gaussian width in z in µm
// Variable vD = wBeta[%D] //diffusion coefficient in µm^2/s
// variable vR = wBeta[%R]//distance from origin in xy plane in µm
// Variable vM = wBeta[%M] //ammount of diffusing substance
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


//copy parameters from wave to local variables for compact expressions
Variable vSigmaXY = wBeta[%SigmaXY] //gaussian width in xy in µm
Variable vSigmaZ = wBeta[%SigmaZ] //gaussian width in z in µm
Variable vD = wBeta[%D] //diffusion coefficient in µm^2/s
variable vR = wBeta[%R]//distance from origin in xy plane in µm
Variable vM = wBeta[%M] //ammount of diffusing substance
Variable vt = wBeta[%tt] //ammount of diffusing substance

make /o/n=(1e4+1) Cx
setscale /i x,(-10*vSigmaXY),(10*vSigmaXY),"µm",Cx
setscale /i y,0,0,"M",Cx
Cx=(vM*exp(((-1)*(x^2))/(4*vD*(vt)+(2*(vSigmaXY^2)))))/(((2*pi)^1.5)*(2*vD*(vt)+(vSigmaXY^2))*sqrt(2*vD*(vt)+(vSigmaZ)^2))
return Cx
end
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// fD3dGaussInit returns the analytical solution to diffusion in 3D with gaussian initial conditions
// initial gaussian is assumed to be symmetric in xy
// the solution assumes z = 0
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
Function D3dGaussT()
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//Data Dictionary
// Variable vSigmaXY = wBeta[%SigmaXY] //gaussian width in xy in µm
// Variable vSigmaZ = wBeta[%SigmaZ] //gaussian width in z in µm
// Variable vD = wBeta[%D] //diffusion coefficient in µm^2/s
// variable vR = wBeta[%R]//distance from origin in xy plane in µm
// Variable vM = wBeta[%M] //ammount of diffusing substance
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
wave wBeta = SetConstants()

//copy parameters from wave to local variables for compact expressions
Variable vSigmaXY = wBeta[%SigmaXY] //gaussian width in xy in µm
Variable vSigmaZ = wBeta[%SigmaZ] //gaussian width in z in µm
Variable vD = wBeta[%D] //diffusion coefficient in µm^2/s
variable vR = wBeta[%R]//distance from origin in xy plane in µm
Variable vM = wBeta[%M] //ammount of diffusing substance

make /o/n=(5e4) Ct
setscale /p x,0,1e-6,"s",Ct
setscale /i y,0,0,"M",Ct
Ct=(vM*exp(((-1)*(vR^2))/(4*vD*(x)+(2*(vSigmaXY^2)))))/(((2*pi)^1.5)*(2*vD*(x)+(vSigmaXY^2))*sqrt(2*vD*(x)+(vSigmaZ)^2))
end


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// GluR_ODE calculates the time derivative of the occupied fraction of receptors
// based on a simplified 2-state model from Destexhe, Mainen, and Sejnowski 1.4.1
//the glutamate dynamics are assumed to be strongly dominated by diffusion
//i.e. glutamate concentration >> receptor concentration so that binding and unbinding
//has negligible effects on free glutamate
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
Function GluR_ODE(pw, tt, yw, dydt)
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//Data Dictionary
// Variable vSigmaXY = pw[%SigmaXY] //gaussian width in xy in µm
// Variable vSigmaZ = pw[%SigmaZ] //gaussian width in z in µm
// Variable vD = pw[%D] //diffusion coefficient in µm^2/s
// variable vR = pw[%R]//distance from origin in xy plane in µm
// Variable vM = pw[%M] //ammount of diffusing substance
// variable vKf = pw[%Kf] // forward reaction rate
// variable vKr = pw[%Kr] // reverse reaction rate
// variable vTau = pw[%Tau] // reverse reaction rate
// Wave pw	//parameter wave
// Variable tt 		// time
// Wave yw	// current y values
// Wave dydt	// y derivatives
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

Wave pw	//parameter wave
Variable tt 		// time
Wave yw	// current y values
Wave dydt	// y derivatives


//copy parameters from wave to local variables for compact expressions
Variable vSigmaXY = pw[%SigmaXY] //gaussian width in xy in µm
Variable vSigmaZ = pw[%SigmaZ] //gaussian width in z in µm
Variable vD = pw[%D] //diffusion coefficient in µm^2/s
variable vR = pw[%R]//distance from origin in xy plane in µm
Variable vM = pw[%M] //ammount of diffusing substance
variable vKf = pw[%Kf]
variable vKr = pw[%Kr]
variable vKdecay = pw[%Kdecay] // reverse reaction rate
variable vGlu
if (vKdecay == 0)
vGlu = ((vM*exp(((-1)*(vR^2))/(4*vD*(tt)+(2*(vSigmaXY^2)))))/(((2*pi)^1.5)*(2*vD*(tt)+(vSigmaXY^2))*sqrt(2*vD*(tt)+(vSigmaZ)^2)))
 else
vGlu = (exp(-vKdecay*tt))*((vM*exp(((-1)*(vR^2))/(4*vD*(tt)+(2*(vSigmaXY^2)))))/(((2*pi)^1.5)*(2*vD*(tt)+(vSigmaXY^2))*sqrt(2*vD*(tt)+(vSigmaZ)^2)))
 endif



dydt[0] = vKf*vGlu *(1 - yw[0]) - vKr*yw[0]
return 0
end

Function SolveGluRt()
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//Data Dictionary
// Variable vSigmaXY = wbeta[%SigmaXY] //gaussian width in xy in µm
// Variable vSigmaZ = wbeta[%SigmaZ] //gaussian width in z in µm
// Variable vD = wbeta[%D] //diffusion coefficient in µm^2/s
// variable vR = wbeta[%R]//distance from origin in xy plane in µm
// Variable vM = wbeta[%M] //ammount of diffusing substance
// variable vKf = wbeta[%Kf] // forward reaction rate
// variable vKr = wbeta[%Kr] // reverse reaction rate
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//if wave wBeta does not exists create it with default values
wave beta = SetConstants()
make /o/n=(beta[%NumPntsReceptorR]) ReceptorR = 0
//setscale /p x,0,1e-6,"s",ReceptorR
setscale /i x,0,beta[%TimeScaleReceptorR],"s",ReceptorR
IntegrateODE GluR_ODE, wBeta, ReceptorR
return 0
end

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// BoundFractionVsDist calculates the fraction of receptor bound to glutamate
// at x = (xmin,...,xmax) for t =(0,tmax)
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

Function BoundFractionVsDist()
wave wBeta = setM()
KillWaves /z ReceptorR_Seq
make /o/n=(wBeta[%nXsteps]) BoundFractionMax
setscale /i x,(wBeta[%Xmin]),(wBeta[%Xmax]),"µm",BoundFractionMax
variable i
for(i=1;i<=wBeta[%nXsteps];i=i+1)	// Initialize variables;continue test
  //wBeta[%R] = wBeta[%Xmin] + (i-1)*(wBeta[%Xmax] - wBeta[%Xmin])/(wBeta[%nXsteps])
  wBeta[%R] = pnt2x(BoundFractionMax,(i-1))
  SolveGluRt()
  BoundFractionMax[i-1] = wavemax(ReceptorR)
  if(!exists("ReceptorR_Seq"))
  duplicate/o ReceptorR, ReceptorR_Seq
  else
  KillWaves /z wTemp
  duplicate/o ReceptorR_seq, wTemp
  KillWaves /z ReceptorR_Seq
  Concatenate {wTemp,ReceptorR},ReceptorR_Seq
  endif

endfor						// Execute body code until continue test is FALSE


end


Function FitGluRxnDiff(w,x) : FitFunc
	Wave w
	Variable x
  wave wBeta
  wave ReceptorR, Rmax
  variable rx, vRmax
	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = m*x+d*x+kf*x+kr*x+xmin*x
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 6
	//CurveFitDialog/ w[0] = M
	//CurveFitDialog/ w[1] = D
	//CurveFitDialog/ w[2] = Kf
	//CurveFitDialog/ w[3] = Kr
	//CurveFitDialog/ w[4] = Xmin
  //CurveFitDialog/ w[5] = Kdecay
  wBeta[%M] = w[0]
  wBeta[%D] = w[1]
  wBeta[%Kf] = w[2]
  wBeta[%Kr] = w[3]
  wBeta[%R] = w[4] + x
  wBeta[%Kdecay] = w[5]
  SolveGluRt()
  rx = wavemax(ReceptorR)
  // if(!exists("Rmax"))
  wBeta[%R] = w[4]
  SolveGluRt()
  // make /n=1 Rmax
  vRmax = wavemax(ReceptorR)
  // endif
	return rx/vRmax
End
