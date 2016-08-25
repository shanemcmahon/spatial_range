#pragma rtGlobals=1
//******************************************************************************
//******************************************************************************
//
//******************************************************************************
//******************************************************************************
Function UserContinue(ctrlName) : ButtonControl
	String ctrlName
	DoWindow/K PauseForUser0
End

//******************************************************************************
//******************************************************************************
//
//******************************************************************************
//******************************************************************************

macro DoUncagingAnalysis()
// helper macro for UncagingAnalysis Function
// the global variables s_wavenames and s_path are not available inside functions, therefore somethings are done first in the macro
	String /g DataWaveList = s_wavenames
  String /g uid =  StringFromList(ItemsInList(s_path, ":")-4, s_path , ":") + "_" + StringFromList(ItemsInList(s_path, ":")-3, s_path , ":") + "_" + StringFromList(ItemsInList(s_path, ":")-2, s_path , ":")+"_ts"+s_path[strlen(s_path)-4,strlen(s_path)-2]
  String /g OutputPathStr = s_path
	uid = uid + s_filename[0,strlen(s_filename)-5]
if(!exists("PointSpacingW"))
	make /o/n=1 PointSpacingW
	PointSpacingW = 100
endif

	UncagingAnalysis(DataWaveList)

endmacro

//******************************************************************************
//******************************************************************************
//
//******************************************************************************
//******************************************************************************

function UserDefineInitialEstimates(ParametersIn,w_coef,UncageTime,y0timeWindow,Amplitude0window, ResponseMaxTime, DelayToResponseStart,FitStop)
Wave ParametersIn, w_coef
Variable UncageTime, &y0timeWindow, &Amplitude0window, &ResponseMaxTime, &DelayToResponseStart, &FitStop
//
//graph average response and get user input for initial estimates
//
//graph
dowindow /k review
display /n=review ParametersIn
SetAxis bottom (UncageTime-y0timeWindow),FitStop
SetDrawEnv xcoord = bottom;SetDrawEnv dash= 3;DelayUpdate
//draw line indicates the location of the uncaging pulse in the aligned average
//because we have aligned the resonse traces at the begining of the fit window, the uncaging pulse occurs at time = UncageTime
DrawLine UncageTime,0,UncageTime,1
ShowInfo/CP=0/W=review
cursor a,$StringFromList(0, tracenamelist("",";",1) ),(UncageTime-y0timeWindow)
//user interaction
//estimate response onset delay
NewPanel/K=2 /n=PauseForUser0 as "Pause for user"; AutoPositionWindow/M=1/R=review
DrawText 21,20,"Adjust cursor A to estimate response start";	DrawText 21,40,"time."
Button button0,pos={80,58},size={92,20},title="Continue"; Button button0,proc=UserContinue
PauseForUser PauseForUser0, review
DelayToResponseStart = xcsr(a)-UncageTime
k1 = xcsr(a)
//estimate peak location and amplitude
NewPanel/K=2 /n=PauseForUser0 as "Pause for user"; AutoPositionWindow/M=1/R=review
DrawText 21,20,"Adjust cursor A to peak response.";//	DrawText 21,40,"Click Continue."
Button button0,pos={80,58},size={92,20},title="Continue"; Button button0,proc=UserContinue
PauseForUser PauseForUser0, review
ResponseMaxTime = xcsr(a)
k0 = mean(ParametersIn,ResponseMaxTime-Amplitude0window,ResponseMaxTime+Amplitude0window) - mean(ParametersIn,0,UncageTime)
k3 = (ResponseMaxTime-k1)*0.33
//estimate decay time
NewPanel/K=2 /n=PauseForUser0 as "Pause for user"; AutoPositionWindow/M=1/R=review
DrawText 21,20,"Adjust cursor A to indicate time at which the";	DrawText 21,40," response has decayed by 90%."
Button button0,pos={80,58},size={92,20},title="Continue"; Button button0,proc=UserContinue
PauseForUser PauseForUser0, review
k2 = (xcsr(a)-ResponseMaxTime)*0.33

k4 = mean(ParametersIn,0,UncageTime)
W_Coef = {k0, DelayToResponseStart, k2, k3, k4, UncageTime}

cursor a,$StringFromList(0, tracenamelist("",";",1) ),(FitStop)


NewPanel/K=2 /n=PauseForUser0 as "Pause for user"; AutoPositionWindow/M=1/R=review
DrawText 21,20,"If necessary, adjust cursor A to edit fit range";	DrawText 21,40,""
Button button0,pos={80,58},size={92,20},title="Continue"; Button button0,proc=UserContinue
PauseForUser PauseForUser0, review
FitStop = xcsr(a)

end

//******************************************************************************
//******************************************************************************
//
//******************************************************************************
//******************************************************************************

function UncagingAnalysis(DataWaveList)
	String DataWaveList
	//name of the uncaging response wave
	String UncagingResponseWaveName;	Prompt UncagingResponseWaveName,"uncaging response wave name",popup,DataWaveList+"Some other wave..."
//name of the uncaging power wave
	String UncagingPowerWaveName; 	Prompt UncagingPowerWaveName,"uncaging power wave name",popup,DataWaveList+"Some other wave..."
//variable that holds the number of uncaging pulses found in UncagingPowerWave
	Variable nUncagingPulses
	//length of the time window used for fitting double exponential
	Variable FitRange = 0.05;	Prompt FitRange,"Size of time window for fitting"
//dummy variable for iteration control
	Variable i = 0, j = 0
		//initial estimate for uncaging response decay time
	Variable DecayTime0 = 0.008;	Prompt DecayTime0,"response decay time initial estimate"
	//initial estimate for uncaging response rise time
	Variable RiseTime0 = 0.001;	Prompt RiseTime0,"response decay time initial estimate"
	Variable DelayToResponseStart = 0
	Variable Amplitude0
	Variable Amplitude0window = 0.001; Prompt Amplitude0window,"Window size for amplitude estimate"
	//time window before uncaging pulse used to estimate y0
	Variable y0timeWindow = 0.01;	Prompt y0timeWindow,"Window size for y0 estimate."
//time of uncaging event for current fit
	Variable UncageTime
	//time of peak amplitude
	Variable ResponseMaxTime
	Variable UserResponse
	//threshold value used while examining UncagingPowerWave to determine whether an algorithmically found peak is a true pulse
	Variable threshold
	Make /O/D w_coef = NaN
	Make /O/D w_sigma = NaN
	//return value from Igor's built-in findlevel function, used for finding uncaging pulses
	Variable V_LevelX = 0
	//initial estimates for fit parameters
	Variable k0, k1, k2, k3, k4
	//variables defining start and stop of the fit window
	Variable FitStart, FitStop
	//Igor environment variable indicating fitting abort condition
	Variable V_AbortCode = 0, V_FitError=0
	// number of points included in fit
	Variable nFitPoints = Inf
	//number of "fake responses" to fit, per stimulus
	Variable NFalseReplicates = 50
	Variable InterFitTime //time between end of a fit and start of the next fit for consecutive fits to "real" stimuli
	Variable InterFalseFitTime //time between fake responses
	Variable ResponseMaxTime0, DelayToResponseStart0 //intial parameter estimates
	Variable Amplitude, TotalLengthUncaging
	// Distance between uncaging points in nm
Wave PointSpacingW
	Variable PointSpacingV = PointSpacingW[0]; prompt PointSpacingV,"Distance between uncaging points"
	//make /o/n=1 PointSpacingW
	wave w_Resampled, wColumnMeans
	make /o /n=1 vPockelsVoltage
	Variable UserSetPar0; prompt UserSetPar0,"Initial parameter estimates",popup,"Interactive;Default;Auto Guess"
	variable v_flag

	Variable LaswerPowerWaveScaling = 0.1; Prompt LaswerPowerWaveScaling, "Laser power wave scaling"


	String WaveDataList =  "AmplitudeW;T0W;DecayTimeW;RiseTimeW;UncageTimeW;OnsetDelayW;y0W;AmplitudeSeW;"

// prompt user to select stimulus and response waves from among recently loaded waves
// user is also provided the option "some other wave..."
// if the user indicates that the stimulus/response wave is not lasted, then prompt again, offering a list of all waves in the current data folder
	Variable nDataWaves = itemsinlist(WaveDataList)
	doPrompt "",UncagingResponseWaveName,UncagingPowerWaveName
	if(!(cmpstr(UncagingResponseWaveName, "Some other wave..." )*cmpstr(UncagingPowerWaveName, "Some other wave..." )))
	Prompt UncagingResponseWaveName,"uncaging response wave name",popup,wavelist("*",";","")
	Prompt UncagingPowerWaveName,"uncaging power wave name",popup,wavelist("*",";","")
	doPrompt "",UncagingResponseWaveName,UncagingPowerWaveName
	endif


// set default value for UserSetPar0
	UserSetPar0 = 3

// promp user for starting parameters
	DoPrompt "",FitRange,DecayTime0,RiseTime0,Amplitude0window,y0timeWindow,UserSetPar0,LaswerPowerWaveScaling,PointSpacingV
	PointSpacingW = PointSpacingV

// set stimulus and response wave references from chosen names
//	wave UncagingResponseWave = $UncagingResponseWaveName
//	wave UncagingPowerWave = $UncagingPowerWaveName


duplicate /o $UncagingResponseWaveName UncagingResponseWave
duplicate /o $UncagingPowerWaveName UncagingPowerWave

	//before we can make waves to put the fit parameters, we need to know how long to make them
	//loop through uncaging events to count the number of uncaging pulses
	//we find the first two pulses manually before entering the do loop
	//as we treat the first pulse slightly different than the rest, treating the first pulses outside of the loop allows us to avoid using an if-then construct inside the loop
	threshold = 0.8*wavemax(UncagingPowerWave)

//calculate average pockels voltage
duplicate /o UncagingPowerWave temp
temp = threshold < UncagingPowerWave
TotalLengthUncaging = sum(temp)
temp = temp * UncagingPowerWave
vPockelsVoltage[0] = sum(temp)/TotalLengthUncaging

	// find rising edge in power trace greater than defined threshold
	findlevel /Q/EDGE=1 /R= (V_LevelX,) UncagingPowerWave, threshold
	//we find the first pulse manually before entering the do loop, so we start the counter at 1
	Make/O/N=1 UncageTimeW
	UncageTimeW = V_LevelX
	i = 1
	do//do1
		// findlevel may occasionaly find false peaks in the noise during the uncaging pulse, but we can reliably find the time where the voltage drops below the threshold
		findlevel /Q/EDGE=2 /R= (V_LevelX,) UncagingPowerWave, threshold
		//the next findlevel searches again for a rising pulse crossing the threshold
		findlevel /Q/EDGE=1 /R= (V_LevelX,) UncagingPowerWave, threshold
		//if no edge is found, V_LevelX will be set to NaN, this condition indicates that we have passed the last uncaging pulse
		if(numtype(V_LevelX))
			break
		endif
		//at this point in the do loop we have succesfully found an uncaging event, so we increment the counter and save the time
		InsertPoints i, 1, UncageTimeW
		UncageTimeW[i] = V_LevelX
		i = i + 1
		//we never have i>100 in experiments so if i>100 something is wrong with the code and we abort
		if(i>100)//if1
			abort("n>100 uncaging pulses found")
		endif//if1
	while(1)//do1
	nUncagingPulses = i
	print "number of uncaging events found: ",nUncagingPulses

	//create waves to store parameters
	make /o/n=(nUncagingPulses) FitStartPointWave
	make /o/n=(nUncagingPulses) FitStartTimeWave
	make /o/n=(nUncagingPulses) FitStopTimeWave
	make /o/n=(nUncagingPulses) AmplitudeW
	make /o/n=(nUncagingPulses) AmplitudeSeW
	make /o/n=(nUncagingPulses) T0W
	make /o/n=(nUncagingPulses) DecayTimeW
	make /o/n=(nUncagingPulses) RiseTimeW
	make /o/n=(nUncagingPulses) y0W
	make /o/n=(nUncagingPulses) OnsetDelayW
	make /o/n=(nUncagingPulses) AmplitudeRestrictedModelWave
	make /o/n=(nUncagingPulses) AmplitudeFromMeanWave
	make /o/n=(nUncagingPulses) AmplitudeRestrictedModelSeWave
	make /o/n=(nUncagingPulses,6) ParametersW2d
	make /o/n=(nUncagingPulses,6) ParameterSeW2d
	//make /o/n=(nUncagingPulses,6) w2dFitUncgRespCoefBootSE

// To get initial estimates of the parameters, the model is first fit to the average response
// to calculate the average response, we need to align the uncaging response relative to the uncaging pulse

//for pre allocating the array to hold the responses, we need to know the number of points, rather than the time
//since the time window is a real number, it's possible that the integral index ranges could have different lengths due to rounding
//if the number of points differs, we use the smallest for analysis
//!note maybe we want to just throw an error if the number of points differs instead
	for(i=0;i < nUncagingPulses;i+=1)	// for1
		FitStart = UncageTimeW[i] - y0timeWindow
		FitStartTimeWave[i] = FitStart
		UncageTime = FitStart + y0timeWindow
		FitStop = FitStart + FitRange
		FitStopTimeWave[i] = FitStop
		nFitPoints = min(nFitPoints,	(x2pnt(UncagingResponseWave, FitStop )-x2pnt(UncagingResponseWave, FitStart )))
		FitStartPointWave[i] = x2pnt(UncagingResponseWave, FitStart )
	endfor		//for1

//copy uncaging responses into 2d wave
make /o/n=(nUncagingPulses, nFitPoints) ResponseWave2d
make /o/n=(nUncagingPulses, nFitPoints) ModelPredictionWave2d
make /o/n=(nFitPoints) ModelPredictionWave
for(i=0;i < nUncagingPulses;i+=1)	// for1
ResponseWave2d[i][] = UncagingResponseWave[FitStartPointWave[i]+q]
endfor		//for1
setscale /p y,0,dimdelta(UncagingResponseWave,0),ResponseWave2d
setscale /p y,0,dimdelta(UncagingResponseWave,0),ModelPredictionWave2d
setscale /p x,0,dimdelta(UncagingResponseWave,0),ModelPredictionWave


//calculate average response
make /o /n=(nUncagingPulses) TempW
TempW = 1
MatrixOp/O ResponseMeanWave=TempW^t x ResponseWave2d
redimension /n=(nFitPoints) ResponseMeanWave
ResponseMeanWave = ResponseMeanWave/nUncagingPulses
setscale /p x,0, dimdelta(UncagingResponseWave,0), ResponseMeanWave

FitStop = FitRange

// ============================================================================
// ============================================================================
// set initial values for parameter estimates
// ============================================================================
// ============================================================================
	if (UserSetPar0 == 1)
	// set user parameters interactively
	// to get user input for initial paramters we call UserDefineInitialEstimates
	// the function has inputs (ParametersIn,w_coef,UncageTime,y0timeWindow,Amplitude0window)
	// //because we have aligned the resonse traces at the begining of the fit window, the uncaging pulse occurs at time = y0timeWindow
	UserDefineInitialEstimates(ResponseMeanWave,w_coef,y0timeWindow,y0timeWindow,Amplitude0window, ResponseMaxTime, DelayToResponseStart,FitStop)
	ResponseMaxTime0 = ResponseMaxTime
	DelayToResponseStart0 = DelayToResponseStart
	V_AbortCode = 0
	V_FitError = 0
	FuncFit/N/Q/H="000001" /NTHR=0 DiffTwoExp2 W_coef  ResponseMeanWave /D
	DelayToResponseStart0 = w_coef[1]
	Amplitude0 = w_coef[0]
	DecayTime0 = w_coef[2]
	RiseTime0 = w_coef[3]
	endif
	if(UserSetPar0 == 2)
	// use parameter values set in dialog
		// ResponseMaxTime = y0timeWindow + 3 * RiseTime0
		DelayToResponseStart = 0
		Amplitude0 = -10
    // duplicate /o /r=(0,y0timeWindow) ResponseMeanWave WTemp
    // wavestats /q WTemp
		w_coef = {-10,DelayToResponseStart,RiseTime0,DecayTime0,mean(ResponseMeanWave,0,y0timeWindow),y0timeWindow}
		// perform fit of average response for display purposes
		FuncFit/N/Q/H="000001" /NTHR=0 DiffTwoExp2 W_coef  ResponseMeanWave /D
		ResponseMaxTime = y0timeWindow + DelayToResponseStart0 + (DecayTime0*RiseTime0)/(DecayTime0-RiseTime0)*ln(DecayTime0/RiseTime0)
		ResponseMaxTime0 = ResponseMaxTime
		DelayToResponseStart = w_coef[1]
		DelayToResponseStart0 = DelayToResponseStart
	endif
	if(UserSetPar0==3)
	// auto guess start parameters
	ResponseMaxTime = y0timeWindow + 3 * RiseTime0
	DelayToResponseStart = 0
	Amplitude0 = -10
	w_coef = {-10,DelayToResponseStart,RiseTime0,DecayTime0,mean(ResponseMeanWave,0,y0timeWindow),y0timeWindow}
	FuncFit/N/Q/H="000001" /NTHR=0 DiffTwoExp2 W_coef  ResponseMeanWave /D
	duplicate /o w_coef wAvgUncageResponseFitCoef
	duplicate /o w_sigma wAvgUncageResponseFitCoefSE
	DelayToResponseStart0 = w_coef[1]
	DelayToResponseStart = DelayToResponseStart0
	Amplitude0 = w_coef[0]
	DecayTime0 = w_coef[2]
	RiseTime0 = w_coef[3]
	ResponseMaxTime = y0timeWindow + w_coef[1] + (w_coef[2]*w_coef[3])/(w_coef[2]-w_coef[3])*ln(w_coef[2]/w_coef[3])
	ResponseMaxTime0 = ResponseMaxTime

	endif


	// ============================================================================
	// ============================================================================
	// do fits
	// ============================================================================
	// ============================================================================

//Make/O/T/N=2 T_Constraints
//T_Constraints[0] = {"K1 > 0","K1 < .01"}

	Prompt UserResponse, "Is this fit good?", popup, "Yes: Save fit;No: Do a Refit; No response: Save zero; Too noisy, save NaN"

	for(i=0;i < nUncagingPulses;i+=1)	// for1
	//the logical control for refitting is handled by performing the fit in a do-while loop with while(UserResponse=2)
	// set user response = 2 to initially enter the loop
	UserResponse = 2

// calculate time window for ith uncaging event
		UncageTime = UncageTimeW[i]
		FitStart = UncageTime - y0timeWindow
		FitStop = FitStart + FitRange

	// set w_coef to initial parameter estimates
		k4 = mean(UncagingResponseWave,FitStart,UncageTime)
		k0 = mean(UncagingResponseWave,(FitStart + ResponseMaxTime0-Amplitude0window),(FitStart + ResponseMaxTime0 + Amplitude0window)) - k4
		k1 = UncageTime + DelayToResponseStart0
		k2 = DecayTime0
		k3 = RiseTime0
    W_Coef = {k0, DelayToResponseStart0, k2, k3, k4, UncageTime}
		V_AbortCode = 0
		V_FitError = 0

// perform initial fit with T0,RiseTime,DecayTime parameters fixed
		// FuncFit/N/Q/H="001101" /NTHR=0 DiffTwoExp2 W_coef  UncagingResponseWave(FitStart, FitStop) /D /C=T_Constraints
		FuncFit/N/Q/H="011101" /NTHR=0 DiffTwoExp2 W_coef  UncagingResponseWave(FitStart, FitStop) /D
		// save the amplitude for restricted model
		AmplitudeRestrictedModelWave[i] = w_coef[0]
		AmplitudeRestrictedModelSeWave[i] = w_sigma[0]
		make /o /n=(nFitPoints) TempW
		TempW = ResponseWave2d[i][p]
		setscale /p x,0, dimdelta(UncagingResponseWave,0), TempW
		duplicate /o TempW wThisResponse
		//save nonparametric estimate of amplitude
		AmplitudeFromMeanWave[i] = mean(TempW,(ResponseMaxTime0-Amplitude0window),(ResponseMaxTime0+Amplitude0window)) - mean(TempW,0,y0timeWindow)
		// save nonparametric estimate of amplitude
		// ResponseMaxTime = y0timeWindow + w_coef[1] + (w_coef[2]*w_coef[3])/(w_coef[2]-w_coef[3])*ln(w_coef[2]/w_coef[3])
		// AmplitudeFromMeanWave[i] = mean(UncagingResponseWave,(FitStart+ResponseMaxTime-Amplitude0window),(FitStart+ResponseMaxTime+Amplitude0window)) - mean(UncagingResponseWave,FitStart,(FitStart+y0timeWindow))

do
// open display window for checking the fit; igor will automatically append the fit to the graph when funcfit is called
dowindow /k review
display /n=review UncagingResponseWave[x2pnt(UncagingResponseWave, FitStart ),x2pnt(UncagingResponseWave, FitStop )]
SetDrawEnv xcoord= bottom;SetDrawEnv dash= 3;DelayUpdate
DrawLine UncageTime,0,UncageTime,1

// perform fit to the full model
// FuncFit/N/Q/H="000001" /NTHR=0 DiffTwoExp2 W_coef  UncagingResponseWave(FitStart, FitStop) /D /C=T_Constraints
FuncFit/N/Q/H="000001" /NTHR=0 DiffTwoExp2 W_coef  UncagingResponseWave(FitStart, FitStop) /D


DoUpdate

// allow user to interact to indicate whether the fit is good

UserResponse = 1 //set default response in the dialog to indicate fit is good
DoPrompt "Goodness of Fit", UserResponse
if(v_flag)
	abort("user canceled")
endif
switch(UserResponse)	// numeric switch
case 1:		// fit is good, nothing to do
	break						// exit from switch
case 2: //user indicated to perform a refit, call UserDefineInitialEstimates
	UserDefineInitialEstimates(UncagingResponseWave,w_coef,UncageTime,y0timeWindow,Amplitude0window, ResponseMaxTime, DelayToResponseStart, FitStop)
	break
case 3: // user indicates no response

	w_coef = 0
	w_sigma = NaN
	break
case 4: //user indicates to save NaN
	w_coef = NaN
	w_sigma = NaN
endswitch

while(UserResponse == 2) //if user indicated to perform a refit, continue loop, else break

// save fit parameters
		AmplitudeW[i] = w_coef[0]
		AmplitudeSeW[i] = w_sigma[0]
		T0W[i] = w_coef[1]
		DecayTimeW[i] = w_coef[2]
		RiseTimeW[i] = w_coef[3]
		y0W[i] = w_coef[4]
		OnsetDelayW[i] = w_coef[1]
		FitStop = FitStart + FitRange
		duplicate /o /r=(FitStart, FitStop) UncagingResponseWave w_t
		w_t = x
		duplicate /o w_t ModelPredictionWave
		ModelPredictionWave = DiffTwoExp2(w_coef, w_t)
		ModelPredictionWave2d[i][] = ModelPredictionWave[q]
		ParametersW2d[i][] = w_coef[q]
		//ParametersW2d[i][] = w_coef[q]
		ParameterSeW2d[i][] = w_sigma[q]





	endfor												// for1; loop through uncaging events, performing fits for each

dowindow/k review



// ============================================================================
// ============================================================================
// perform fits to response trace over periods where no stimulus is given
// ============================================================================
// ============================================================================

make /o /n=((nUncagingPulses-1)*NFalseReplicates) NullAmpRestricedW
	make /o /n=((nUncagingPulses-1)*NFalseReplicates) NullAmplitudeFullModelWave
make /o /n=((nUncagingPulses-1)*NFalseReplicates) NullAmplitudeFromMeanWave

for(i=0;i < (nUncagingPulses-1);i+=1)	// for1
 // there is no stimulus between uncaging pulses
 // loop over uncaging pulses and perform fits during the time between uncaging events
InterFitTime = FitStartTimeWave[i+1] - FitStopTimeWave[i]
InterFalseFitTime = InterFitTime/NFalseReplicates
for(j=0;j < NFalseReplicates;j+=1)	// for2
// during each inter-stimulus, perform NFalseReplicates fits

// set initial parameters and fit window
FitStart = FitStopTimeWave[i] + j*InterFalseFitTime
UncageTime = FitStart + y0timeWindow
FitStop = FitStart + FitRange
k4 = mean(UncagingResponseWave,FitStart,UncageTime)
k0 = mean(UncagingResponseWave,(FitStart + ResponseMaxTime0 - Amplitude0window),(FitStart + ResponseMaxTime0 - Amplitude0window)) - k4
Amplitude = k0
k1 = UncageTime + DelayToResponseStart
k2 = DecayTime0
k3 = RiseTime0
W_Coef = {k0, DelayToResponseStart0, k2, k3, k4,UncageTime}
V_AbortCode = 0
V_FitError = 0

// fit restricted model
FuncFit/N/Q/W=2/H="011101" /NTHR=0 DiffTwoExp2 W_coef  UncagingResponseWave(FitStart, FitStop) /D
// FuncFit/N/Q/W=2/H="001101" /NTHR=0 DiffTwoExp2 W_coef  UncagingResponseWave(FitStart, FitStop) /D /C=T_Constraints
// ResponseMaxTime = y0timeWindow + w_coef[1] + (w_coef[2]*w_coef[3])/(w_coef[2]-w_coef[3])*ln(w_coef[2]/w_coef[3])
// save nonparametric estimate of amplitude
// NullAmplitudeFromMeanWave[i*NFalseReplicates+j] = mean(UncagingResponseWave,(FitStart+ResponseMaxTime-Amplitude0window),(FitStart+ResponseMaxTime+Amplitude0window)) - mean(UncagingResponseWave,FitStart,(FitStart+y0timeWindow))
NullAmplitudeFromMeanWave[i*NFalseReplicates+j] = Amplitude
// save amplitude from restricted model
NullAmpRestricedW[i*NFalseReplicates+j] = w_coef[0]
// perform full model fit
// FuncFit/N/Q/W=2/H="000001" /NTHR=0 DiffTwoExp2 W_coef  UncagingResponseWave(FitStart, FitStop) /D /C=T_Constraints
FuncFit/N/Q/W=2/H="000001" /NTHR=0 DiffTwoExp2 W_coef  UncagingResponseWave(FitStart, FitStop) /D
// duplicate /o/r=[0,4] w_coef w_coef2
// duplicate /o /r =(FitStart, FitStop) UncagingResponseWave Wtemp
// duplicate /o wTemp wTemp2
// wTemp2 = x
// concatenate /o {wTemp2,wTemp}, ModelFrame
// optimize /q /x=w_coef /R=w_coef /m={0,1} DiffTwoExpSM,w_coef
// w_coef[0,4]=w_coef2[p]
// save amplitude from full model fit
NullAmplitudeFullModelWave[i*NFalseReplicates+j] = w_coef[0]
// NullAmplitudeFullModelWave[i*NFalseReplicates+j] = w_coef[0]

// save nonparametric amplitude estimate
// NullAmplitudeFromMeanWave[i*NFalseReplicates+j] = Amplitude

endfor//for2
endfor//for1

end


//******************************************************************************
//******************************************************************************
//
//******************************************************************************
//******************************************************************************

Function DiffTwoExp2(w,t) : FitFunc
	Wave w
	Variable t



	IF    ((t-(w[1]+w[5])) < 0 )
		return w[4]
	ELSE
		//Difference in exponential model from Schutter, Erik De. Computational modeling methods for neuroscientists. The MIT Press, 2009. Chapter 6
		 variable Gsyn, Fnorm, TPeak,t0,td,tr,y0,UncageTime
		Gsyn = w[0]
		t0 = w[1]
		td = w[2]
		tr = w[3]
		y0 = w[4]
		UncageTime = w[5]
		TPeak = (UncageTime+t0) + (td*tr)/(td-tr)*ln(td/tr)
		fnorm = 1/(-exp(-(Tpeak-(UncageTime+t0))/tr) + exp(-(Tpeak-(UncageTime+t0))/td))
		return (y0 + Gsyn*fnorm*(exp(-(t-(UncageTime+t0))/td)-exp(-(t-(UncageTime+t0))/tr)))
    // return (w[4]+w[0]*(-exp(-(t-(w[1]+w[5]))/(w[3]))+exp(-(t-(w[1]+w[5]))/(w[2])))/(-exp(-((w[1])+((w[2])*(w[3])/((w[2])-(w[3])))*ln((w[2])/(w[3]))-(w[1]))/(w[3]))+exp(-((w[1])+((w[2])*(w[3])/((w[2])-(w[3])))*ln((w[2])/(w[3]))-(w[1]))/(w[2]))))
	ENDIF
End

//******************************************************************************
//******************************************************************************
//
//******************************************************************************
//******************************************************************************

Function DiffTwoExpSM(par,xIn)
	Wave par,xIn
	wave ModelFrame, w_coef


	variable Gsyn, Fnorm, TPeak,t0,td,tr,y0,UncageTime
	Gsyn = par[0]
	t0 = par[1]
	td = par[2]
	tr = par[3]
	y0 = par[4]
	UncageTime = par[5]

	//if(dimsize(ModelFrame,0 < dimsize(ModelFrame,1)))
	if(dimsize(ModelFrame,0) < dimsize(ModelFrame,1))
	matrixop /o xTemp = ModelFrame^t
	ModelFrame = xTemp
	endif

	make /o /n=(dimsize(ModelFrame,0)) Yhat
	make /o /n=(dimsize(ModelFrame,0)) yObs
	make /o /n=(dimsize(ModelFrame,0)) tIn
	tIn = ModelFrame[p][0]
	yObs = ModelFrame[p][1]

	make /o /n=(numpnts(tIn)) TMaskBaseline
	make /o /n=(numpnts(tIn)) TMaskNotBaseline
	TMaskBaseline = tIn < (t0+UncageTime)
	TMaskNotBaseline = !TMaskBaseline
		TPeak = (UncageTime+t0) + (td*tr)/(td-tr)*ln(td/tr)
		fnorm = 1/(-exp(-(Tpeak-(UncageTime+t0))/tr) + exp(-(Tpeak-(UncageTime+t0))/td))

matrixop/o yhat = (y0)*TMaskBaseline + TMaskNotBaseline*(y0 + Gsyn*fnorm*(exp(-(tIn-(UncageTime+t0))/td)-exp(-(tIn-(UncageTime+t0))/tr)))
// duplicate /o yobs errs
matrixop /o errs = yobs - Yhat
matrixop /o sse = errs.errs
// errs = errs*errs
variable sse_out = sse[0]
//print sse_out,par
return sse_out
		//return (y0 + Gsyn*fnorm*(exp(-(t-(UncageTime+t0))/td)-exp(-(t-(UncageTime+t0))/tr)))
    // return (w[4]+w[0]*(-exp(-(t-(w[1]+w[5]))/(w[3]))+exp(-(t-(w[1]+w[5]))/(w[2])))/(-exp(-((w[1])+((w[2])*(w[3])/((w[2])-(w[3])))*ln((w[2])/(w[3]))-(w[1]))/(w[3]))+exp(-((w[1])+((w[2])*(w[3])/((w[2])-(w[3])))*ln((w[2])/(w[3]))-(w[1]))/(w[2]))))

End

//******************************************************************************
//******************************************************************************
//
//******************************************************************************
//******************************************************************************


macro do_save_results()
if(!exists("rw_uid"))
make /o /n=0 /t rw_uid
endif
InsertPoints numpnts(rw_uid), 1, rw_uid
rw_uid[numpnts(rw_uid)] = uid
save_results()
endmacro

//******************************************************************************
//******************************************************************************
//
//******************************************************************************
//******************************************************************************

function save_results()
wave rw2d_response, UncagingResponseWave,UncageTimeW
wave FitStartTimeWave, FitStopTimeWave, AmplitudeW, AmplitudeSeW, T0W
wave DecayTimeW, RiseTimeW, y0W, OnsetDelayW, AmplitudeRestrictedModelWave
wave AmplitudeFromMeanWave, ResponseWave2d, ModelPredictionWave2d, rw_uid, NullAmpRestricedW
wave NullAmplitudeFullModelWave, NullAmplitudeFromMeanWave, rwPockelsVoltage
wave vPockelsVoltage
variable n_results
string /g OutputPathStr
string /g uid
newpath /o OutputDir, OutputPathStr
SavePICT/O/E=-5/B=288 /p=OutputDir /win=SummaryFig as (uid +".png")
//SavePICT/O/E=-5/B=72 /p=OutputDir /win=SummaryFig as "SummaryFiglr.png"

if(!waveexists(rw2d_response))

make /n=(numpnts(UncagingResponseWave),8) rw2d_response
make /n=(numpnts(UncageTimeW),8) rw2d_uncage_time
make /n=(numpnts(FitStartTimeWave),8) rw2d_fit_start_time
make /n=(numpnts(FitStopTimeWave),8) rw2d_fit_stop_time
make /n=(numpnts(AmplitudeW),8) rw2d_fit_amplitude
make /n=(numpnts(AmplitudeSeW),8) rw2d_fit_amplitude_se
make /n=(numpnts(T0W),8) rw2d_fit_t0
make /n=(numpnts(DecayTimeW),8) rw2d_fit_decay_time
make /n=(numpnts(RiseTimeW),8) rw2d_fit_rise_time
make /n=(numpnts(y0W),8) rw2d_fit_y0
make /n=(numpnts(OnsetDelayW),8) rw2d_fit_onset_delay
make /n=(numpnts(AmplitudeRestrictedModelWave),8) rw2d_amplitude_0
make /n=(numpnts(AmplitudeFromMeanWave),8) rw2d_amplitude_0_np
// make /n=(numpnts(),8)


duplicate ResponseWave2d rw3d_uncaging_response
redimension /n=(-1,-1,8) rw3d_uncaging_response
duplicate ModelPredictionWave2d rw3d_fits
redimension /n=(-1,-1,8) rw3d_fits
duplicate NullAmpRestricedW W2dNrAmplitude0
redimension /n=(-1,8) W2dNrAmplitude0
duplicate NullAmplitudeFullModelWave W2dNrAmplitude1
redimension /n=(-1,8) W2dNrAmplitude1
duplicate NullAmplitudeFromMeanWave W2dNrAmplitude2
redimension /n=(-1,8) W2dNrAmplitude2
make /o /n=0 WAmplitudeCorrelation
make /o /n=0 WAmplitudeNrCorrelation
make /o /n=0 rwPockelsVoltage
endif

n_results = numpnts(rw_uid)-1

if(n_results >  dimsize( rw2d_response, 1)*3/4)
Redimension /N=(-1, 2*n_results) rw2d_response
Redimension /N=(-1, 2*n_results) rw2d_uncage_time
Redimension /N=(-1, 2*n_results) rw2d_fit_start_time
Redimension /N=(-1, 2*n_results) rw2d_fit_stop_time
Redimension /N=(-1, 2*n_results) rw2d_fit_amplitude
Redimension /N=(-1, 2*n_results) rw2d_fit_amplitude_se
Redimension /N=(-1, 2*n_results) rw2d_fit_t0
Redimension /N=(-1, 2*n_results) rw2d_fit_decay_time
Redimension /N=(-1, 2*n_results) rw2d_fit_rise_time
Redimension /N=(-1, 2*n_results) rw2d_fit_y0
Redimension /N=(-1, 2*n_results) rw2d_fit_onset_delay
Redimension /N=(-1, 2*n_results) rw2d_amplitude_0
Redimension /N=(-1, 2*n_results) rw2d_amplitude_0_np
// Redimension /N=(-1, 2*n_results)

Redimension /N=(-1,-1, 2*n_results) rw3d_uncaging_response
Redimension /N=(-1,-1, 2*n_results) rw3d_fits
redimension /n=(-1,2*n_results) W2dNrAmplitude0
redimension /n=(-1,2*n_results) W2dNrAmplitude1
redimension /n=(-1,2*n_results) W2dNrAmplitude2


endif


rw2d_response[][n_results] = UncagingResponseWave[p]
rw2d_uncage_time[][n_results] = UncageTimeW[p]
rw2d_fit_start_time[][n_results] = FitStartTimeWave[p]
rw2d_fit_stop_time[][n_results] = FitStopTimeWave[p]
rw2d_fit_amplitude[][n_results] = AmplitudeW[p]
rw2d_fit_amplitude_se[][n_results] = AmplitudeSeW[p]
rw2d_fit_t0[][n_results] = T0W[p]
rw2d_fit_decay_time[][n_results] = DecayTimeW[p]
rw2d_fit_rise_time[][n_results] = RiseTimeW[p]
rw2d_fit_y0[][n_results] = y0W[p]
rw2d_fit_onset_delay[][n_results] = OnsetDelayW[p]
rw2d_amplitude_0[][n_results] = AmplitudeRestrictedModelWave[p]
rw2d_amplitude_0_np[][n_results] = AmplitudeFromMeanWave[p]
W2dNrAmplitude0[][n_results] = NullAmpRestricedW[p]
W2dNrAmplitude1[][n_results] = NullAmplitudeFullModelWave[p]
W2dNrAmplitude2[][n_results] = NullAmplitudeFromMeanWave[p]
// [][n_results] = [p]

rw3d_uncaging_response[][][n_results]=ResponseWave2d[p][q]
rw3d_fits[][][n_results]=ModelPredictionWave2d[p][q]

InsertPoints numpnts(WAmplitudeCorrelation), 1, WAmplitudeCorrelation
InsertPoints numpnts(WAmplitudeNrCorrelation), 1, WAmplitudeNrCorrelation
InsertPoints numpnts(rwPockelsVoltage), 1, rwPockelsVoltage
rwPockelsVoltage[n_results] = vPockelsVoltage[0]
WAmplitudeCorrelation[n_results] = StatsCorrelation(AmplitudeW,AmplitudeFromMeanWave)
WAmplitudeNrCorrelation[n_results] = StatsCorrelation(NullAmplitudeFullModelWave,NullAmplitudeFromMeanWave)
end

//******************************************************************************
//******************************************************************************
//
//******************************************************************************
//******************************************************************************

macro Clean_Up()
	dowindow/k graph5
	dowindow/k graph4
	dowindow/k graph3
	dowindow/k graph2
	dowindow/k graph1
	dowindow/k graph0
  dowindow /k review
	//killwindow layout0
	Kill_Input_Waves()
kill_wave_list("ACH_1;ACH_3;")
kill_wave_list("UncagingResponseWave;UncagingPowerWave;")
kill_wave_list("fit_w2d_responses;w_t;w2d_fake_pars;fit_w_response_out;w_bs_amp0;w_bs_amp0alt;w_bs_amp0_Hist;w_bs_amp0alt_Hist;")
kill_wave_list("ResponseWave2d;ModelPredictionWave2d;ModelPredictionWave;TempW;ResponseMeanWave;T_Constraints;fit_w_avg_response;")
kill_wave_list("RiseTimeW;y0W;OnsetDelayW;AmplitudeRestrictedModelWave;AmplitudeFromMeanWave;AmplitudeRestrictedModelSeWave;")
kill_wave_list("FitStartTimeWave;FitStopTimeWave;AmplitudeW;AmplitudeSeW;T0W;DecayTimeW;")
kill_wave_list("w_response_out;w_power_out;w_coef;W_sigma;UncageTimeW;FitStartPointWave;")
kill_wave_list("ACH_1;ACH_3;UncageTimeW;w_refs;w_stim1;w_stim2;w_stim3;w_stim4;w_stim5;w_response_out;w_power_out;ResponseWave2d;w2d_stim;TempW;ResponseMeanWave;w_avg_power;")
	killstrings /a/z
endmacro

//******************************************************************************
//******************************************************************************
//
//******************************************************************************
//******************************************************************************

macro Kill_Input_Waves()
	//
	//Kill_Input_Waves kills the waves listed in the DataWaveList String
	//if DataWaveList does not exist, then the function will attempt to kill the waves listed in s_wavenames
	//s_wavenames is auto generated by igor whenever a data wave is loaded
	//this function would normally be invoked by the user after DoUncagingAnalysis which would kill only the input waves but leave the output from the analysis in place
	//

	if(!exists("DataWaveList"))//if1
		if(!exists("s_wavenames"))//if2
			return
		endif//endif2
		String DataWaveList = s_wavenames
	endif//endif1
	Kill_Wave_List(DataWaveList)
endmacro

//******************************************************************************
//******************************************************************************
//
//******************************************************************************
//******************************************************************************


function RemoveSpineData()
variable i_
wave rw2d_fit_amplitude,rw2d_amplitude_0,rw2d_amplitude_0_np,rw2d_fit_amplitude_se
wave rw2d_fit_decay_time,rw2d_fit_onset_delay,rw2d_fit_rise_time,rw2d_fit_start_time,rw2d_fit_stop_time
wave rw2d_fit_t0,rw2d_fit_y0,rw2d_response,rw2d_uncage_time,rw3d_fits,rw3d_uncaging_response
wave rwPockelsVoltage,rw_uid

prompt i_,"Point number"
doprompt "Enter value",i_

DeletePoints/M=1 i_,1, rw2d_fit_amplitude
DeletePoints/M=1 i_,1, rw2d_amplitude_0
DeletePoints/M=1 i_,1, rw2d_amplitude_0_np
DeletePoints/M=1 i_,1, rw2d_fit_amplitude_se
DeletePoints/M=1 i_,1, rw2d_fit_decay_time
DeletePoints/M=1 i_,1, rw2d_fit_onset_delay
DeletePoints/M=1 i_,1, rw2d_fit_rise_time
DeletePoints/M=1 i_,1, rw2d_fit_start_time
DeletePoints/M=1 i_,1, rw2d_fit_stop_time
DeletePoints/M=1 i_,1, rw2d_fit_t0
DeletePoints/M=1 i_,1, rw2d_fit_y0
DeletePoints/M=1 i_,1, rw2d_response
DeletePoints/M=1 i_,1, rw2d_uncage_time
DeletePoints/M=2 i_,1, rw3d_fits
DeletePoints/M=2 i_,1, rw3d_uncaging_response
DeletePoints i_,1, rwPockelsVoltage
DeletePoints i_,1, rw_uid
end

//******************************************************************************
//******************************************************************************
//
//******************************************************************************
//******************************************************************************

function SetResponseNaN()
wave AmplitudeW,AmplitudeSeW,T0W,DecayTimeW,RiseTimeW,y0W,OnsetDelayW
wave AmplitudeRestrictedModelWave,AmplitudeFromMeanWave,AmplitudeRestrictedModelSeWave,ResponseWave2d,ResponseWave2d,ModelPredictionWave2d
variable i_
prompt i_,"Point number"
doprompt "Enter value",i_
//ShowInfo/CP=0
//cursor a,$StringFromList(0, tracenamelist("",";",1) ),0
AmplitudeW[i_]=NaN;AmplitudeSeW[i_]=NaN;T0W[i_]=NaN;DecayTimeW[i_]=NaN;RiseTimeW[i_]=NaN;y0W[i_]=NaN;
OnsetDelayW[i_]=NaN;AmplitudeRestrictedModelWave[i_]=NaN;AmplitudeFromMeanWave[i_]=NaN;AmplitudeRestrictedModelSeWave[i_]=NaN;
ResponseWave2d[i_][]=NaN
ModelPredictionWave2d[i_][]=NaN
end

//******************************************************************************
//******************************************************************************
//
//******************************************************************************
//******************************************************************************

function RemoveResponse()
wave AmplitudeW,AmplitudeSeW,T0W,DecayTimeW,RiseTimeW,y0W,OnsetDelayW
wave AmplitudeRestrictedModelWave,AmplitudeFromMeanWave,AmplitudeRestrictedModelSeWave,ResponseWave2d,ResponseWave2d,ModelPredictionWave2d
variable i_
prompt i_,"Point number"
doprompt "Enter value",i_
DeletePoints i_,1, AmplitudeW,AmplitudeSeW,T0W,DecayTimeW,RiseTimeW,y0W;DelayUpdate
DeletePoints i_,1, OnsetDelayW,AmplitudeRestrictedModelWave,AmplitudeFromMeanWave,AmplitudeRestrictedModelSeWave
DeletePoints i_,1, ResponseWave2d
DeletePoints i_,1, ModelPredictionWave2d
end

//******************************************************************************
//******************************************************************************
//
//******************************************************************************
//******************************************************************************

function InsertResponse()
wave AmplitudeW,AmplitudeSeW,T0W,DecayTimeW,RiseTimeW,y0W,OnsetDelayW
wave AmplitudeRestrictedModelWave,AmplitudeFromMeanWave,AmplitudeRestrictedModelSeWave,ResponseWave2d,ResponseWave2d,ModelPredictionWave2d
Insertpoints 0,1, AmplitudeW,AmplitudeSeW,T0W,DecayTimeW,RiseTimeW,y0W;DelayUpdate
Insertpoints 0,1, OnsetDelayW,AmplitudeRestrictedModelWave,AmplitudeFromMeanWave,AmplitudeRestrictedModelSeWave
Insertpoints 0,1, ResponseWave2d
Insertpoints 0,1, ModelPredictionWave2d

AmplitudeW[0]=NaN;AmplitudeSeW[0]=NaN;T0W[0]=NaN;DecayTimeW[0]=NaN;RiseTimeW[0]=NaN;y0W[0]=NaN;
OnsetDelayW[0]=NaN;AmplitudeRestrictedModelWave[0]=NaN;AmplitudeFromMeanWave[0]=NaN;AmplitudeRestrictedModelSeWave[0]=NaN;
ResponseWave2d[0][]=NaN
ModelPredictionWave2d[0][]=NaN

end

//******************************************************************************
//******************************************************************************
//
//******************************************************************************
//******************************************************************************

macro DoMakeFigures(UncageSpacingV)
	//MakeFigures()
	variable UncageSpacingV = PointSpacingW[0]
	variable vNumStim

	variable i
	String cmd

	//vPockelsVoltage = 20*vPockelsVoltage
	dowindow/k FullResponse
	display /n=FullResponse UncagingResponseWave
	Label bottom "time \\u#2 (s)"
	Label left "I\\u#2 (pA)"
	wavestats /	q UncagingResponseWave
	SetAxis left v_min,(v_max+v_sdev)
	duplicate /o UncageTimeW wUncageIndicator
	wUncageIndicator = (v_max+v_sdev)
	appendtograph wUncageIndicator vs UncageTimeW
	ModifyGraph mode(wUncageIndicator)=3,marker(wUncageIndicator)=2;DelayUpdate
	ModifyGraph rgb(wUncageIndicator)=(0,0,0)
	ModifyGraph width={Aspect,7}
	ModifyGraph width=430,height=225
	ModifyGraph margin(bottom)=40
	ModifyGraph margin(left)=40
	ModifyGraph margin(top)=20
	ModifyGraph margin(right)=20

	vNumStim = dimsize(ResponseWave2d,0)
	duplicate /o UncageTimeW wUncageIndicatorTime
	wUncageIndicatorTime = UncageTimeW[0]-FitStartTimeWave[0]
	i = 0
	do
	sprintf cmd, "dowindow /k response%s", num2str(i)
	print cmd
	Execute cmd
	sprintf cmd, "Display /n= response%s ResponseWave2d[%s][*]", num2str(i), num2str(i)
	print cmd
	Execute cmd
		Label bottom "\\u#2time (ms)"
		Label left "I\\u#2 (pA)"
		appendtograph ModelPredictionWave2d[i][*]
		ModifyGraph rgb(ModelPredictionWave2d)=(0,0,0)
		appendtograph wUncageIndicator vs wUncageIndicatorTime
		ModifyGraph mode(wUncageIndicator)=3,marker(wUncageIndicator)=2;DelayUpdate
		ModifyGraph rgb(wUncageIndicator)=(0,0,0)
		SetAxis left v_min,(v_max+v_sdev)
		ModifyGraph width=185,height=80
		ModifyGraph margin=50
		ModifyGraph margin(top)=0,margin(right)=0
		ModifyGraph lblMargin=10

		i += 1
	while (i < vNumStim)






	dowindow/k FitAmplitude
	display /n=FitAmplitude AmplitudeW
	setscale /p x,((numpnts(AmplitudeW)-1)*UncageSpacingV),-UncageSpacingV,"nm",AmplitudeW
	Label bottom "distance \\u#2 (nm)"
	ModifyGraph mode=3,marker=19;DelayUpdate
//	ErrorBars/T=0/L=0.7 AmplitudeW Y,wave=(AmplitudeSeW,AmplitudeSeW)
	Label left "I\\u#2 (pA)"
	SetAxis left *,max(wavemax(AmplitudeW),0)
	Sort NullAmplitudeFullModelWave NullAmplitudeFullModelWave
	setscale /i x,0,1,NullAmplitudeFullModelWave
	SetDrawEnv ycoord= left;SetDrawEnv dash= 3;DelayUpdate
	DrawLine 0,(NullAmplitudeFullModelWave(0.01)),1,(NullAmplitudeFullModelWave(0.01))
	SetDrawEnv ycoord= left;SetDrawEnv dash= 3;DelayUpdate
	DrawLine 0,(NullAmplitudeFullModelWave(0.05)),1,(NullAmplitudeFullModelWave(0.05))
	SetDrawEnv ycoord= left;SetDrawEnv dash= 3;DelayUpdate
	DrawLine 0,(NullAmplitudeFullModelWave(0.5)),1,(NullAmplitudeFullModelWave(0.5))
	ModifyGraph width=185,height=80
	ModifyGraph margin=50
	ModifyGraph margin(top)=0,margin(right)=0
	ModifyGraph lblMargin=10
	//AutoPositionWindow/M=1/R=graph1

	dowindow/k UncageResponses
	plotrows(ResponseWave2d)
	dowindow /c UncageResponses
	ModifyGraph width=185,height=80
	ModifyGraph margin=50
	ModifyGraph margin(top)=0,margin(right)=0
	ModifyGraph lblMargin=10

	//AutoPositionWindow/M=0/R=graph4 graph5

	dowindow/k UncageFits
	plotrows(ModelPredictionWave2d)
	dowindow /c UncageFits
	ModifyGraph width=185,height=80
	ModifyGraph margin=50
	ModifyGraph margin(top)=0,margin(right)=0
	ModifyGraph lblMargin=10
	//AutoPositionWindow/M=0/R=graph5 graph6
	duplicate /o VPockelsVoltage VLaserPower
	if(!exists("LaserPower"))
		make LaserPower
	endif
	setscale /p x, 0, 0.05, LaserPower
	VLaserPower = LaserPower(VPockelsVoltage)

	make/o/t/n=(3,2) SummaryFigTableData
	summaryFigTableData[][0] = {"Pockels (mV)","Laser (mW)","ID"}
	summaryFigTableData[][1] = {"","",""}
	summaryFigTableData[0][1] = {num2str(vPockelsVoltage[0]),num2str(VLaserPower[0]),uid}
	dowindow /k SummaryFigTable
	edit /n=SummaryFigTable SummaryFigTableData
	modifytable autosize = {0,1,-1,0,0}
	ModifyTable width(Point)=0
	modifytable alignment=0

MakeLayout()

endmacro


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
macro MakeLayout(nResponsePanelCols)
Variable nResponsePanelCols = 2
Variable i,j
String cmd

dowindow /k SummaryFig
newlayout /n=SummaryFig
ModifyGraph /w=FullResponse width=(245*nResponsePanelCols-60)
appendtolayout FullResponse
modifylayout left(FullResponse)=0
modifylayout top(FullResponse)=0

duplicate /o AmplitudeW UncagePosition
UncagePosition = x


i=0

do
print i
j = 0
do
//print j
//print (i*nResponsePanelCols+j)
sprintf cmd, "appendtolayout Response%s",num2str(i*nResponsePanelCols+j);	Execute cmd
sprintf cmd, "modifylayout left(Response%s)=%s",num2str(i*nResponsePanelCols+j),num2str(245*j);	Execute cmd
sprintf cmd, "modifylayout top(Response%s)=%s",num2str(i*nResponsePanelCols+j),num2str(285 + i*135);	Execute cmd
sprintf cmd, "TextBox/w=SummaryFig/C/N=text%s/F=0 /A=LT (num2str(UncagePosition[%s]) + \"nm\")",num2str(i*nResponsePanelCols+j),num2str(i*nResponsePanelCols+j);	Execute cmd
sprintf cmd, "modifylayout left(text%s)=%s",num2str(i*nResponsePanelCols+j),num2str(245*j);	Execute cmd
sprintf cmd, "modifylayout top(text%s)=%s",num2str(i*nResponsePanelCols+j),num2str(285 + i*135);	Execute cmd
j += 1
while(j<nResponsePanelCols)


if(i*nResponsePanelCols + j + 1 >= dimsize(ModelPredictionWave2d,0))
break
endif

i += 1
while(i < dimsize(ModelPredictionWave2d,0))

appendtolayout UncageResponses
modifylayout left(UncageResponses)=0
modifylayout top(UncageResponses)=420+(i)*135
appendtolayout UncageFits
modifylayout left(UncageFits)=245
modifylayout top(UncageFits)=420+(i)*135

appendtolayout FitAmplitude
modifylayout left(FitAmplitude)=0
modifylayout top(FitAmplitude)=535+(i)*135

appendtolayout SummaryFigTable
modifylayout left(SummaryFigTable)=245
modifylayout top(SummaryFigTable)=535+(i)*135
modifylayout frame = 0
endmacro



menu "macros"
	"DoUncagingAnalysis/1"
	"DoMakeFigures/2"
	"Do_Save_Results/3"
	"Clean_Up/4"
	"Kill_Input_Waves"
	"SetResponseNan/5"
	"RemoveResponse/6"
	"InsertResponse/7"
	"RemoveSpineData"
	"MakeLayout"
end
