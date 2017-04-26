#pragma rtGlobals=1
#include ":utility_functions"
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
	newdatafolder /o UAdata
	newdatafolder /o FitResults

	//generate spine id. If s_filename exists, then it is used as the id, otherwise the current data folder is used
	make /t /o/n=1 uid
	if(exists("s_filename"))
		uid = s_filename[0,strlen(s_filename)-5]
	else
		uid = getdatafolder(0)
	endif

	//if s_path does not exists, then call pathinfo home, this sets s_path to the home directory
	if(!exists("s_path"))
		pathinfo home
	endif
	String /g OutputPathStr = s_path
	make /o /t /n=1 wOutputPathStr = s_path
	//provide wave containers and default values for some user specefied parameters
	//if the waves already exist then they are not overwritten, this allows parameters
	//set by user to be preserved between function calls

	if(!waveexists(:UAdata:PointSpacingW))
		make /o/n=1 :UAdata:PointSpacingW
		:UAdata:PointSpacingW = -100
	endif

	if(!waveexists(:UAdata:X0W))
		make /o/n=1 :UAdata:X0W
		:UAdata:X0W = -1000
	endif

	//delay between uncaging and response (T0) at distance = 0 nm from spine
	if(!waveexists(:UAdata:T00W))
		make /o/n=1 :UAdata:T00W
		:UAdata:T00W = 6.6E-4
	endif
	//change in T0 with distance in s per nm
	if(!waveexists(:UAdata:DeltaT0W))
		make /o/n=1 :UAdata:DeltaT0W
		:UAdata:DeltaT0W = 1.1E-6
	endif
	UncagingAnalysis()
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

function UncagingAnalysis()
	//******************************************************************************
	//******************************************************************************
	//declare variables, waves, and strings
	//name of the uncaging response wave
	String UncagingResponseWaveName;	Prompt UncagingResponseWaveName,"uncaging response wave name",popup,wavelist("*",";","")+"Some other wave..."
	//name of the uncaging power wave
	String UncagingPowerWaveName; 	Prompt UncagingPowerWaveName,"uncaging power wave name",popup,wavelist("*",";","")+"Some other wave..."
	//variable that holds the number of uncaging pulses found in UncagingPowerWave
	Variable nUncagingPulses
	//length of the time window used for fitting double exponential
	Variable FitRange = 0.05;	Prompt FitRange,"Size of time window for fitting(s)"
	//dummy variable for iteration control
	Variable i = 0, j = 0
	//initial estimate for uncaging response decay time
	Variable DecayTime0 = 0.008;	Prompt DecayTime0,"response decay time initial estimate (s)"
	//initial estimate for uncaging response rise time
	Variable RiseTime0 = 0.001;	Prompt RiseTime0,"response rise time initial estimate (s)"
	Variable DelayToResponseStart = 0
	//initial estimate for response amplitude
	Variable Amplitude0
	Variable Amplitude0window = 0.001; Prompt Amplitude0window,"Window size for amplitude estimate(s)"
	//time window before uncaging pulse used to estimate y0
	Variable y0timeWindow = 0.01;	Prompt y0timeWindow,"Window size for y0 estimate(s)"
	//time of uncaging event for current fit
	Variable UncageTime
	//time of peak amplitude
	Variable ResponseMaxTime
	//dummy variable to hold user response to prompts
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
	Variable ResponseMaxTime0, DelayToResponseStart0 //intial parameter estimates
	Variable TotalLengthUncaging
	// Distance between uncaging points in nm, stored in a wave for persistance
	Wave PointSpacingW = :UAdata:PointSpacingW
	Variable PointSpacingV = PointSpacingW[0]; prompt PointSpacingV,"Distance between uncaging points (nm)"
	Wave T00W = :UAdata:T00W
	Variable T00V = T00W[0]; prompt T00V,"T00"
	Wave DeltaT0W = :UAdata:DeltaT0W
	Variable DeltaT0V = DeltaT0W[0]; prompt DeltaT0V,"DeltaT0"
	//wColumnMeans is the output of a user defined function (ColumnMeans), to make igor happy we need to tell it that the wave exists
	make /o wColumnMeans

	//bit string specifying which model parameters are fixed/free. Used as funcfit /H=CoefIsFixed
	String CoefIsFixed = "011101"; prompt CoefIsFixed, "Bit string specifying parameters to fix in the model." //see Igor help on funcfit /H for details
	make /o /n=1 vPockelsVoltage
	Variable UserSetPar0; prompt UserSetPar0,"Initial parameter estimates",popup,"Interactive;Default;Auto Guess"
	variable v_flag
	//CurrentWaveNames is a wave containing the names of the waves used in the current instance (response, power, etc.)
	make /o/t /n=(1,10) :UAdata:CurrentWaveNames
	wave /t CurrentWaveNames = :UAdata:CurrentWaveNames
	CurrentWaveNames = ""
	wave /t uid
	string suid = uid[0]; prompt suid,"Terminal Unique ID"
	string HelpString


	// ============================================================================
	// ============================================================================
	// query user for inputs
	// ============================================================================
	// ============================================================================



	String WaveDataList =  "AmplitudeW;T0W;DecayTimeW;RiseTimeW;UncageTimeW;OnsetDelayW;y0W;AmplitudeSeW;"

	// prompt user to select stimulus and response waves from among recently loaded waves
	// user is also provided the option "some other wave..."
	// if the user indicates that the stimulus/response wave is not lasted, then prompt again, offering a list of all waves in the current data folder
	Variable nDataWaves = itemsinlist(WaveDataList)
	HelpString = "UncagingResponseWaveName: the name of the wave containing the uncaging response. The units are assumed to be seconds and amps.\n\n UncagingPowerWaveName: the name of the wave containing the pockels voltage. The units are assumed to be in seconds and volts"
	doPrompt /help=HelpString "",UncagingResponseWaveName,UncagingPowerWaveName

	CurrentWaveNames[0][0] = UncagingResponseWaveName; SetDimLabel 1,0,response,CurrentWaveNames
	CurrentWaveNames[0][1] = UncagingPowerWaveName; SetDimLabel 1,1,power,CurrentWaveNames


	// set default value for UserSetPar0
	UserSetPar0 = 3

	// promp user for starting parameters
	HelpString = "Size of time window: Time window in seconds for response fittings.\n\n"
	HelpString = HelpString + "Response Decay Time: initial estimate for the response decay time, in seconds\n\n"
	HelpString = HelpString + "Response rise time: initial estimate for response rise time, in seconds\n\n"
	HelpString = HelpString + "Window size for amplitude estimate: The response is averaged over a time window with the specified width. The center of the time window is calculated from the rise and decay times in auto mode. In interactive mode, the center is specefied interactively."
	HelpString = HelpString + "Window size for y0 estimate: The response is averaged over the specified time interval to estimate the initial amplitude. In auto mode the start of the window is set to 10 ms before the uncaging pulse. In interactive mode, it can be specefied manually."
	HelpString = HelpString + "Initial paramter estimates: \n 	Defualt: uses default parameters.\n Auto Guess: The program attempts to find suitable starting parameters automatically. \n Interactive: Estimate starting parameters by interacting with a graph.\n\n"
	HelpString = HelpString + "Distance between uncaging points: The distance between uncaging points. The sign indicates direction of motion. A negative value indicates that consective uncages pulses approach the spine. A positive value indicates that succesive pulses are moving away from the spine. The sign of the value determines which response traces are used for estimating initial parameters.\n\n"
	HelpString = HelpString + "T00: The parameter T0 supplied to the Difference in exponential fitting function is the time between the uncaging pulse and the rising phase of the response. T00 sets T0 for the response closest to the spine"
	HelpString = HelpString + "DeltaT0: T0 is assumed to increase linearly with spine distance according to T0 = T00 + DeltaT0*Distance. DeltaT0 has units of ms/nm"

	DoPrompt /Help="Refer to manual for parameter descriptions." "",FitRange,DecayTime0,RiseTime0,Amplitude0window,y0timeWindow,UserSetPar0,PointSpacingV,T00V,DeltaT0V

	HelpString = 	"CoefIsFixed: bit string specifying which model parameters are fixed/free. Used as funcfit /H=CoefIsFixed\n\n"
	HelpString = HelpString + "suid: a name to identify the response. The spine id is used to name output files."
	DoPrompt /Help=HelpString "", CoefIsFixed, suid
	uid[0] = suid
	PointSpacingW = PointSpacingV
	T00W = T00V
	DeltaT0W = DeltaT0V

	// set stimulus and response wave references from chosen names

	wave UncagingResponseWave = $UncagingResponseWaveName
	wave UncagingPowerWave = $UncagingPowerWaveName


	// ============================================================================
	// ============================================================================
	// find uncaging pulses and allocate variables
	// ============================================================================
	// ============================================================================
	//before we can make waves to put the fit parameters, we need to know how long to make them
	//loop through uncaging events to count the number of uncaging pulses
	//we find the first two pulses manually before entering the do loop
	//as we treat the first pulse slightly different than the rest, treating the first pulses outside of the loop allows us to avoid using an if-then construct inside the loop
	threshold = 0.8*wavemax(UncagingPowerWave)

	//calculate average pockels voltage
	duplicate /o UncagingPowerWave temp
	// create a mask to define when uncaging power is greater than threshold
	temp = threshold < UncagingPowerWave
	// calculate average when mask is true
	TotalLengthUncaging = sum(temp)
	temp = temp * UncagingPowerWave
	vPockelsVoltage[0] = sum(temp)/TotalLengthUncaging*20

	// find rising edge in power trace greater than defined threshold
	findlevel /Q/EDGE=1 /R= (V_LevelX,) UncagingPowerWave, threshold
	//we find the first pulse manually before entering the do loop, so we start the counter at 1
	Make/O/N=1 :UAdata:UncageTimeW
	wave UncageTimeW = :UAdata:UncageTimeW
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
		//we never have i>100 in experiments so if i>100 something is wrong with the code or input data and we abort
		if(i>100)//if1
			abort("n>100 uncaging pulses found")
		endif//if1
	while(1)//do1
	nUncagingPulses = i
	print "number of uncaging events found: ",nUncagingPulses

	//create waves to store parameters
	make /o/n=(nUncagingPulses) :UAdata:FitStartPointWave; wave FitStartPointWave = :UAdata:FitStartPointWave
	make /o/n=(nUncagingPulses) :UAdata:FitStartTimeWave; wave FitStartTimeWave = :UAdata:FitStartTimeWave
	make /o/n=(nUncagingPulses) :UAdata:FitStopTimeWave; wave FitStopTimeWave = :UAdata:FitStopTimeWave
	make /o/n=(nUncagingPulses) :FitResults:AmplitudeW; wave AmplitudeW = :FitResults:AmplitudeW
	make /o/n=(nUncagingPulses) :FitResults:AmplitudeSeW; wave AmplitudeSeW = :FitResults:AmplitudeSeW
	make /o/n=(nUncagingPulses) :FitResults:T0W; wave T0W = :FitResults:T0W
	make /o/n=(nUncagingPulses) :FitResults:DecayTimeW; wave DecayTimeW = :FitResults:DecayTimeW
	make /o/n=(nUncagingPulses) :FitResults:RiseTimeW; wave RiseTimeW = :FitResults:RiseTimeW
	make /o/n=(nUncagingPulses) :FitResults:y0W; wave y0W = :FitResults:y0W
	make /o/n=(nUncagingPulses) :FitResults:OnsetDelayW; wave OnsetDelayW = :FitResults:OnsetDelayW


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
	make /o/n=(nUncagingPulses, nFitPoints) :FitResults:ResponseWave2d; wave ResponseWave2d = :FitResults:ResponseWave2d
	make /o/n=(nUncagingPulses, nFitPoints) :FitResults:ModelPredictionWave2d; wave ModelPredictionWave2d = :FitResults:ModelPredictionWave2d
	make /o/n=(nFitPoints) ModelPredictionWave
	for (i=0;i < nUncagingPulses;i+=1)	// for1
		ResponseWave2d[i][] = UncagingResponseWave[FitStartPointWave[i]+q]
		Duplicate /O /r =[i][] ResponseWave2d  $(":UAdata:posNo_"+num2str(i))
		redimension /n=(nfitpoints) $(":UAdata:posNo_"+num2str(i))
		setscale /p x,0,dimdelta(UncagingResponseWave,0), $(":UAdata:posNo_"+num2str(i))
	endfor		//for1
	setscale /p y,0,dimdelta(UncagingResponseWave,0),ResponseWave2d
	setscale /p y,0,dimdelta(UncagingResponseWave,0),ModelPredictionWave2d
	setscale /p x,0,dimdelta(UncagingResponseWave,0),ModelPredictionWave

	// ============================================================================
	// ============================================================================
	// calculate average over large responses used to set fit starting parameters
	// ============================================================================
	// ============================================================================





	//if user supplies a wave then use it for estimating starting parameters
	//if user supplied wave name is empty then average the first three uncaging responses
	//if user supplies a wave name but the reference does not exists, then abort

	make /o/n=(3,nFitPoints) :UAdata:BigResponseW
	wave BigResponseW = :UAdata:BigResponseW
	BigResponseW[][] = ResponseWave2d[nUncagingPulses-p-1][q]
	if(PointSpacingV> 0)
		BigResponseW[][] = ResponseWave2d[p][q]
	endif
	ColumnMeans(BigResponseW)
	duplicate /o wColumnMeans, ResponseMeanWave


	String UserAvgResponseWaveName = ""; prompt UserAvgResponseWaveName, "Name of user specified average wave name (optional)",popup,wavelist("*",";","")
	UserAvgResponseWaveName = "ResponseMeanWave"
	DoPrompt /HELP="Specify the wave name used to estimate initial parameters. \n If ResponseMeanWave (the default value) is specified, then the program estimates initial parameters using an average of 3 responses, else the program uses the named wave." "", UserAvgResponseWaveName
	// doprompt "", UserAvgResponseWaveName
	duplicate /o $UserAvgResponseWaveName temp
	duplicate /o temp ResponseMeanWave
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
		//because we have aligned the resonse traces at the begining of the fit window, the uncaging pulse occurs at time = y0timeWindow
		UserDefineInitialEstimates(ResponseMeanWave,w_coef,y0timeWindow,y0timeWindow,Amplitude0window, ResponseMaxTime, DelayToResponseStart,FitStop)
		ResponseMaxTime0 = ResponseMaxTime
		DelayToResponseStart0 = DelayToResponseStart
		V_AbortCode = 0
		V_FitError = 0
		Display /N=look ResponseMeanWave
		ModifyGraph rgb(ResponseMeanWave)=(0,0,0)
		FuncFit/N/Q/H="000001" /NTHR=0 DiffTwoExp2 W_coef  ResponseMeanWave /D

		NewPanel/K=2 /n=PauseForUser0 as "Pause for user"; AutoPositionWindow/M=1/R=look
		Button button0,pos={80,58},size={92,20},title="Continue"; Button button0,proc=UserContinue
		PauseForUser PauseForUser0, look
		dowindow /k look

		//save parameters
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
	if(UserSetPar0 == 2)
		// use parameter values set in dialog
		DelayToResponseStart = 0
		Amplitude0 = -10
		w_coef = {-10,DelayToResponseStart,RiseTime0,DecayTime0,mean(ResponseMeanWave,0,y0timeWindow),y0timeWindow}
		// perform fit of average response for display purposes
		Display /N=look ResponseMeanWave
		ModifyGraph rgb(ResponseMeanWave)=(0,0,0)
		FuncFit/N/Q/H="000001" /NTHR=0 DiffTwoExp2 W_coef  ResponseMeanWave /D
		//review fit
		NewPanel/K=2 /n=PauseForUser0 as "Pause for user"; AutoPositionWindow/M=1/R=look
		Button button0,pos={80,58},size={92,20},title="Continue"; Button button0,proc=UserContinue
		PauseForUser PauseForUser0, look
		dowindow /k look
		//save parameters
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
	if(UserSetPar0==3)
		// auto guess start parameters
		//fit averaged responses
		ResponseMaxTime = y0timeWindow + 3 * RiseTime0
		DelayToResponseStart = 0
		Amplitude0 = -10
		w_coef = {-10,DelayToResponseStart,RiseTime0,DecayTime0,mean(ResponseMeanWave,0,y0timeWindow),y0timeWindow}
		Display /N=look ResponseMeanWave
		ModifyGraph rgb(ResponseMeanWave)=(0,0,0)
		FuncFit/N/Q/H="000001" /NTHR=0 DiffTwoExp2 W_coef  ResponseMeanWave /D

		//review fit
		NewPanel/K=2 /n=PauseForUser0 as "Pause for user"; AutoPositionWindow/M=1/R=look
		Button button0,pos={80,58},size={92,20},title="Continue"; Button button0,proc=UserContinue
		PauseForUser PauseForUser0, look
		dowindow /k look

		//save parameters
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
		k1 = T00V + DeltaT0V*(i*PointSpacingV)
		DelayToResponseStart0 = k1
		if(PointSpacingV < 0)
			k1 = T00V + DeltaT0V*(abs(PointSpacingV)*(nUncagingPulses-1) + i*PointSpacingV )
			DelayToResponseStart0 = k1
		endif
		if(!DeltaT0V)
			k1 = DelayToResponseStart0
		endif

		k2 = DecayTime0; k3 = RiseTime0
		W_Coef = {k0, k1, k2, k3, k4, UncageTime}
		V_AbortCode = 0;	V_FitError = 0



		make /o /n=(nFitPoints) :UAdata:wThisResponse
		wave wThisResponse = :UAdata:wThisResponse
		wThisResponse = ResponseWave2d[i][p]
		setscale /p x,0, dimdelta(UncagingResponseWave,0), wThisResponse


		do
			// open display window for checking the fit; igor will automatically append the fit to the graph when funcfit is called
			dowindow /k review
			display /n=review UncagingResponseWave[x2pnt(UncagingResponseWave, FitStart ),x2pnt(UncagingResponseWave, FitStop )]


			SetDrawEnv xcoord= bottom;SetDrawEnv dash= 3;DelayUpdate
			DrawLine UncageTime,0,UncageTime,1

			FuncFit/N/Q/H=CoefIsFixed /NTHR=0 DiffTwoExp2 W_coef  UncagingResponseWave(FitStart, FitStop) /D
			DoUpdate

			//******************************************************************************
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

		//******************************************************************************
		// save fit parameters and estimates
		AmplitudeW[i] = w_coef[0]; AmplitudeSeW[i] = w_sigma[0]; T0W[i] = w_coef[1]
		DecayTimeW[i] = w_coef[2]; RiseTimeW[i] = w_coef[3];	y0W[i] = w_coef[4]
		OnsetDelayW[i] = w_coef[1];	FitStop = FitStart + FitRange
		//calculate model estimates
		duplicate /o /r=(FitStart, FitStop) UncagingResponseWave :UAdata:w_t
		wave w_t = :UAdata:w_t
		w_t = x
		duplicate /o w_t ModelPredictionWave
		ModelPredictionWave = DiffTwoExp2(w_coef, w_t)
		//save estimates, parameter matrix
		ModelPredictionWave2d[i][] = ModelPredictionWave[q]
		// ParametersW2d[i][] = w_coef[q];ParameterSeW2d[i][] = w_sigma[q]
	endfor												// for1; loop through uncaging events, performing fits for each

	dowindow/k review

	//******************************************************************************
	//******************************************************************************
	// If point spacing > 0 then uncaging position is moving away from Spine
	// then reverse order of data
	make /o /n=1 UncagingOrderIsReversed
	UncagingOrderIsReversed = 0
	if(PointSpacingV > 0)
		duplicate /o AmplitudeW, TempW;		  AmplitudeW = TempW[(nUncagingPulses-1)-p]
		duplicate /o AmplitudeSeW, TempW;		AmplitudeSeW = TempW[(nUncagingPulses-1)-p]
		duplicate /o T0W, TempW;		        T0W = TempW[(nUncagingPulses-1)-p]
		duplicate /o DecayTimeW, TempW;		  DecayTimeW = TempW[(nUncagingPulses-1)-p]
		duplicate /o RiseTimeW, TempW;		  RiseTimeW = TempW[(nUncagingPulses-1)-p]
		duplicate /o y0W, TempW;		        y0W = TempW[(nUncagingPulses-1)-p]
		duplicate /o OnsetDelayW, TempW;		OnsetDelayW = TempW[(nUncagingPulses-1)-p]

		duplicate /o ModelPredictionWave2d TempW
		ModelPredictionWave2d[][] = TempW[(nUncagingPulses-1)-p][q]
		duplicate /o ResponseWave2d TempW
		ResponseWave2d[][] = TempW[(nUncagingPulses-1)-p][q]

		UncagingOrderIsReversed = 1
	endif

	killwaves /z :UAdata:fit_responseMeanWave, :UAdata:ResponseMeanWave, :UAdata:wAvgUncageResponseFitCoef
	killwaves /z :UAdata:wAvgUncageResponseFitCoefSE, w_coef, w_sigma
	MoveWave fit_responseMeanWave, :UAdata:
	MoveWave ResponseMeanWave, :UAdata:
	MoveWave wAvgUncageResponseFitCoef, :UAdata:
	MoveWave wAvgUncageResponseFitCoefSE, :UAdata:
	killwaves /z ModelPredictionWave, temp,wColumnMeans
	//killwaves /z fit_responseMeanWave
	String fldrSav0= GetDataFolder(1)
	SetDataFolder :FitResults:

	make /o /n=(DimSize(ModelPredictionWave2d, 1 ),DimSize(ModelPredictionWave2d, 0 )) ModelPredictionWave2dT
	ModelPredictionWave2dT = ModelPredictionWave2d[q][p]
	setscale /p y,0,DimDelta(ModelPredictionWave2d, 1 ),ModelPredictionWave2dT
	KillDataFolder /Z FitWaves
	SplitWave/z/O/SDIM=1/N=w_fit/DDF=FitWaves ModelPredictionWave2dT

	make /o /n=(DimSize(ResponseWave2d, 1 ),DimSize(ResponseWave2d, 0 )) ResponseWave2dT
	ResponseWave2dT = ResponseWave2d[q][p]
	setscale /p y,0,DimDelta(ResponseWave2d, 1 ),ResponseWave2dT
	KillDataFolder /Z ResponseWaves
	SplitWave/z/O/SDIM=1/N=w_Response/DDF=ResponseWaves ResponseWave2dT


	SetDataFolder fldrSav0


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


//******************************************************************************
//******************************************************************************
//
//******************************************************************************
//******************************************************************************


macro DoSaveResults()
	SaveResults()
endmacro

//******************************************************************************
//******************************************************************************
//
//******************************************************************************
//******************************************************************************

function SaveResults()
	wave /t CurrentWaveNames = :UAdata:CurrentWaveNames
	wave ResponseW2d
	wave UncageTimeW = :UAdata:UncageTimeW
	wave AmplitudeW = :FitResults:AmplitudeW
	wave AmplitudeSeW = :FitResults:AmplitudeSeW
	wave T0W = :FitResults:T0W
	wave DecayTimeW = :FitResults:DecayTimeW
	wave RiseTimeW = :FitResults:RiseTimeW
	wave y0W = :FitResults:y0W
	wave OnsetDelayW = :FitResults:OnsetDelayW
	wave ResponseWave2d =:FitResults:ResponseWave2d
	wave ModelPredictionWave2d =	:FitResults:ModelPredictionWave2d
	wave PockelsVoltageW
	wave vPockelsVoltage
	wave UncagingResponseWave = $(currentwavenames[0][%response])
	wave /t uid
	duplicate uid, :FitResults:uid

	if(!exists("root:UidW"))
		make /o /n=0 /t root:UidW
	endif
	wave /t UidW = root:UidW
	InsertPoints numpnts(UidW), 1, root:UidW

	root:UidW[numpnts(UidW)] = uid[0]


	variable nResults

	if(!waveexists(root:ResponseW2d))

		make /n=(numpnts(UncagingResponseWave),8) root:ResponseW2d
		make /n=(numpnts(UncageTimeW),8) root:UncageTimeW2d
		make /n=(numpnts(AmplitudeW),8) root:AmplitudeW2d
		make /n=(numpnts(AmplitudeSeW),8) root:AmplitudeSeW2d
		make /n=(numpnts(T0W),8) root:T0W2d
		make /n=(numpnts(DecayTimeW),8) root:DecayTimeW2d
		make /n=(numpnts(RiseTimeW),8) root:RiseTimeW2d
		make /n=(numpnts(y0W),8) root:y0W2d
		make /n=(numpnts(OnsetDelayW),8) root:OnsetDelayW2d



		duplicate ResponseWave2d root:UncagingResponseW3d
		redimension /n=(-1,-1,8) root:UncagingResponseW3d
		duplicate ModelPredictionWave2d root:ModelPredictionW3d
		redimension /n=(-1,-1,8) root:ModelPredictionW3d
		make /o /n=0 root:PockelsVoltageW
	endif
	wave ResponseW2d = root:ResponseW2d
	wave UncageTimeW2d = root:UncageTimeW2d
	wave AmplitudeW2d = root:AmplitudeW2d
	wave AmplitudeSeW2d = root:AmplitudeSeW2d
	wave T0W2d = root:T0W2d
	wave DecayTimeW2d = root:DecayTimeW2d
	wave RiseTimeW2d = root:RiseTimeW2d
	wave y0W2d = root:y0W2d
	wave OnsetDelayW2d = root:OnsetDelayW2d
	wave UncagingResponseW3d = root:UncagingResponseW3d
	wave ModelPredictionW3d = root:ModelPredictionW3d
	wave PockelsVoltageW = root:PockelsVoltageW

	nResults = numpnts(root:UidW)-1

	if(nResults >  dimsize( root:ResponseW2d, 1)*3/4)
		Redimension /N=(-1, 2*nResults) root:ResponseW2d
		Redimension /N=(-1, 2*nResults) root:UncageTimeW2d
		Redimension /N=(-1, 2*nResults) root:AmplitudeW2d
		Redimension /N=(-1, 2*nResults) root:AmplitudeSeW2d
		Redimension /N=(-1, 2*nResults) root:T0W2d
		Redimension /N=(-1, 2*nResults) root:DecayTimeW2d
		Redimension /N=(-1, 2*nResults) root:RiseTimeW2d
		Redimension /N=(-1, 2*nResults) root:y0W2d
		Redimension /N=(-1, 2*nResults) root:OnsetDelayW2d

		Redimension /N=(-1,-1, 2*nResults) root:UncagingResponseW3d
		Redimension /N=(-1,-1, 2*nResults) root:ModelPredictionW3d
	endif



	root:ResponseW2d[][nResults] = UncagingResponseWave[p]
	root:UncageTimeW2d[][nResults] = UncageTimeW[p]
	root:AmplitudeW2d[][nResults] = AmplitudeW[p]
	root:AmplitudeSeW2d[][nResults] = AmplitudeSeW[p]
	root:T0W2d[][nResults] = T0W[p]
	root:DecayTimeW2d[][nResults] = DecayTimeW[p]
	root:RiseTimeW2d[][nResults] = RiseTimeW[p]
	root:y0W2d[][nResults] = y0W[p]
	root:OnsetDelayW2d[][nResults] = OnsetDelayW[p]
	SetDimLabel 1,nResults,$uid[0],AmplitudeW2d
	root:UncagingResponseW3d[][][nResults]=ResponseWave2d[p][q]
	root:ModelPredictionW3d[][][nResults]=ModelPredictionWave2d[p][q]
	InsertPoints numpnts(PockelsVoltageW), 1, root:PockelsVoltageW
	root:PockelsVoltageW[nResults] = vPockelsVoltage[0]

end

//******************************************************************************
//******************************************************************************
//
//******************************************************************************
//******************************************************************************

macro Clean_Up()
	//killwaves UncagingResponseWave
endmacro

//******************************************************************************
//******************************************************************************
//
//******************************************************************************
//******************************************************************************


//******************************************************************************
//******************************************************************************
//
//******************************************************************************
//******************************************************************************


function RemoveSpineData()
	variable i_
	wave AmplitudeW2d,AmpRestrictedModelW2d,AmpFromMeanW2d,AmplitudeSeW2d
	wave DecayTimeW2d,OnsetDelayW2d,RiseTimeW2d,UncageTimeW2dFitStartTimeW2d,FitStopTimeW2d
	wave T0W2d,y0W2d,ResponseW2d,UncageTimeW2d,ModelPredictionW3d,UncagingResponseW3d
	wave PockelsVoltageW,UidW

	prompt i_,"Point number"
	doprompt "Enter value",i_

	DeletePoints/M=1 i_,1, AmplitudeW2d
	DeletePoints/M=1 i_,1, AmplitudeSeW2d
	DeletePoints/M=1 i_,1, DecayTimeW2d
	DeletePoints/M=1 i_,1, OnsetDelayW2d
	DeletePoints/M=1 i_,1, RiseTimeW2d
	DeletePoints/M=1 i_,1, UncageTimeW2d
	DeletePoints/M=1 i_,1, T0W2d
	DeletePoints/M=1 i_,1, y0W2d
	DeletePoints/M=1 i_,1, ResponseW2d
	DeletePoints/M=2 i_,1, ModelPredictionW3d
	DeletePoints/M=2 i_,1, UncagingResponseW3d
	DeletePoints i_,1, PockelsVoltageW
	DeletePoints i_,1, UidW
end

//******************************************************************************
//******************************************************************************
//
//******************************************************************************
//******************************************************************************

function SetResponseNaN()
	wave AmplitudeW,AmplitudeSeW,T0W,DecayTimeW,RiseTimeW,y0W,OnsetDelayW
	wave ResponseWave2d,ResponseWave2d,ModelPredictionWave2d
	variable i_
	prompt i_,"Point number"
	doprompt "Enter value",i_

	AmplitudeW[i_]=NaN;AmplitudeSeW[i_]=NaN;T0W[i_]=NaN;DecayTimeW[i_]=NaN;RiseTimeW[i_]=NaN;y0W[i_]=NaN;
	OnsetDelayW[i_]=NaN;
	ModelPredictionWave2d[i_][]=NaN

end

//******************************************************************************
//******************************************************************************
//
//******************************************************************************
//******************************************************************************

function RemoveResponse()
	wave AmplitudeW,AmplitudeSeW,T0W,DecayTimeW,RiseTimeW,y0W,OnsetDelayW
	wave ResponseWave2d,ResponseWave2d,ModelPredictionWave2d
	variable i_
	prompt i_,"Point number"
	doprompt "Enter value",i_
	DeletePoints i_,1, AmplitudeW,AmplitudeSeW,T0W,DecayTimeW,RiseTimeW,y0W;DelayUpdate
	DeletePoints i_,1, OnsetDelayW
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
	wave ResponseWave2d,ResponseWave2d,ModelPredictionWave2d
	Insertpoints 0,1, AmplitudeW,AmplitudeSeW,T0W,DecayTimeW,RiseTimeW,y0W;DelayUpdate
	Insertpoints 0,1, OnsetDelayW
	Insertpoints 0,1, ResponseWave2d
	Insertpoints 0,1, ModelPredictionWave2d


	AmplitudeW[0]=NaN;AmplitudeSeW[0]=NaN;T0W[0]=NaN;DecayTimeW[0]=NaN;RiseTimeW[0]=NaN;y0W[0]=NaN;
	OnsetDelayW[0]=NaN;
	ResponseWave2d[0][]=NaN
	ModelPredictionWave2d[0][]=NaN
end

//******************************************************************************
//******************************************************************************
//
//******************************************************************************
//******************************************************************************



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
function MakeLayout()
	string CellRefImageWave
	wave /t wCellRefImageWaveName
	Prompt CellRefImageWave,"cell reference image wave name",popup,(wavelist("*",";","dims:2")+wavelist("*",";","dims:3"))

	doprompt "",CellRefImageWave
	wCellRefImageWaveName[0] = CellRefImageWave
end

macro DoMakeLayout(nResponsePanelCols,LayoutWidth,LayoutHeight)
	Variable nResponsePanelCols = nResponsePanelColsW[0]
	Variable LayoutWidth = 1350
	Variable LayoutHeight = 1200

	Variable i,j
	String cmd
	make /O/T wCellRefImageWaveName
	nResponsePanelColsW[0] = nResponsePanelCols

	MakeLayout()
	dowindow /k CellRefImage
	NewImage/K=0/n=CellRefImage $wCellRefImageWaveName[0]
	ModifyGraph /w=CellRefImage tick=3,noLabel=2

	dowindow /k SummaryFig
	newlayout /n=SummaryFig
	LayoutPageAction size(-1)=(LayoutWidth, LayoutHeight), margins(-1)=(18, 18, 18, 18)

	//if(!wintype("SummaryFig"))
	//	newlayout /n=SummaryFig
	//else
	//layoutpageaction /w = SummaryFig appendpage
	//endif


	ModifyGraph /w=FullResponse width=(245*nResponsePanelCols-60)
	appendtolayout FullResponse
	modifylayout left(FullResponse)=0
	modifylayout top(FullResponse)=0

	duplicate /o :FitResults:AmplitudeW UncagePosition
	UncagePosition = x

	ModifyLayout units=0
	i=0
	string ResponseName, TextBoxName
	string layout_info = layoutinfo("SummaryFig","layout");string LayoutSize = stringfromlist(1,stringfromlist(3,layout_info),":")
	variable layoutright = str2num(stringfromlist(2,LayoutSize,",")), variable layoutbottom = str2num(stringfromlist(3,LayoutSize,","))
	do
		print i
		j = 0
		do
			ResponseName = ("Response" + num2str(i*nResponsePanelCols+j))
			TextBoxName = ("TextBox"+num2str(i*nResponsePanelCols+j))
			AppendToLayout $ResponseName
			//AppendToLayout Response1
			// sprintf cmd, "appendtolayout Response%s",num2str(i*nResponsePanelCols+j);	Execute cmd
			ModifyLayout left($ResponseName)=(245*j)
			// sprintf cmd, "modifylayout left(Response%s)=%s",num2str(i*nResponsePanelCols+j),num2str(245*j);	Execute cmd
			ModifyLayout top($ResponseName)=((285 + i*135))
			// sprintf cmd, "modifylayout top(Response%s)=%s",num2str(i*nResponsePanelCols+j),num2str(285 + i*135);	Execute cmd
			TextBox /w=SummaryFig/C/N=$TextBoxName/F=0/A=LT (num2str(UncagePosition[(i*nResponsePanelCols+j)]) + "nm")
			// sprintf cmd, "TextBox/w=SummaryFig/C/N=text%s/F=0 /A=LT (num2str(UncagePosition[%s]) + \"nm\")",num2str(i*nResponsePanelCols+j),num2str(i*nResponsePanelCols+j);	Execute cmd
			ModifyLayout  left($TextBoxName)=(245*j)
			// sprintf cmd, "modifylayout left(text%s)=%s",num2str(i*nResponsePanelCols+j),num2str(245*j);	Execute cmd
			ModifyLayout  top($TextBoxName)=((285 + i*135))
			// sprintf cmd, "modifylayout top(text%s)=%s",num2str(i*nResponsePanelCols+j),num2str(285 + i*135);	Execute cmd
			j += 1
			if(i*nResponsePanelCols + j >= dimsize(:FitResults:ModelPredictionWave2d,0))
				break
			endif
		while(j<nResponsePanelCols)
		if(i*nResponsePanelCols + j >= dimsize(:FitResults:ModelPredictionWave2d,0))
			break
		endif
		i += 1
	while(i < dimsize(:FitResults:ModelPredictionWave2d,0))



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

	appendtolayout CellRefImage
	modifylayout left(CellRefImage)=(245*nResponsePanelCols)
	//	layoutpageaction /w = SummaryFig size=( ((420+(i)*135)+min(512,420+(i)*135)),numberbykey("Top",layoutinfo("SummaryFig","SummaryFigTable")) + numberbykey("Height",layoutinfo("SummaryFig","SummaryFigTable")))


endmacro
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

macro DoMakeFigures(UncageSpacingV,OutputFormat,sTemp)
	variable UncageSpacingV = :UAdata:PointSpacingW[0]
	string OutputFormat = "pdf"
	prompt OutputFormat,"Output Format",popup "png;pdf;svg;"
	string sTemp
	Prompt sTemp, "Output Path", popup (pathlist("*",";","")+"other;")
	duplicate /o $(:UAdata:currentwavenames[0][%response]) UncagingResponseWave

	UncageSpacingV = abs(UncageSpacingV)
	variable vNumStim

	variable i
	String cmd

	//vPockelsVoltage = 20*vPockelsVoltage
	dowindow/k FullResponse
	display /n=FullResponse UncagingResponseWave
	Label bottom "time \\u#2 (s)"; Label left "I\\u#2 (pA)"
	wavestats /	q UncagingResponseWave
	SetAxis left v_min,(v_max+v_sdev)
	duplicate /o :UAdata:UncageTimeW :UAdata:wUncageIndicator
	:UAdata:wUncageIndicator = (v_max+v_sdev)
	appendtograph :UAdata:wUncageIndicator vs :UAdata:UncageTimeW
	ModifyGraph mode(wUncageIndicator)=3,marker(wUncageIndicator)=2;DelayUpdate
	ModifyGraph rgb(wUncageIndicator)=(0,0,0)//, width={Aspect,7}
	ModifyGraph width=430,height=225
	ModifyGraph margin(bottom)=40, margin(left)=40, margin(top)=20, margin(right)=20
	// ModifyGraph margin(left)=40
	// ModifyGraph margin(top)=20
	// ModifyGraph margin(right)=20

	vNumStim = dimsize(:FitResults:ResponseWave2d,0)
	duplicate /o :UAdata:UncageTimeW :UAdata:wUncageIndicatorTime
	:UAdata:wUncageIndicatorTime = UncageTimeW[0]-FitStartTimeWave[0]
	i = 0
	do
		sprintf cmd, "dowindow /k response%s", num2str(i)
		//	print cmd
		Execute cmd
		sprintf cmd, "Display /n= response%s :FitResults:ResponseWave2d[%s][*]", num2str(i), num2str(i)
		//	print cmd
		Execute cmd
		Label bottom "\\u#2time (ms)"
		Label left "I\\u#2 (pA)"
		appendtograph :FitResults:ModelPredictionWave2d[i][*]
		ModifyGraph rgb(ModelPredictionWave2d)=(0,0,0)
		SetDrawEnv xcoord=bottom, ycoord=rel
		DrawText 10E-3,0.1,"*"
		// appendtograph :UAdata:wUncageIndicator vs :UAdata:wUncageIndicatorTime
		// ModifyGraph mode(wUncageIndicator)=3,marker(wUncageIndicator)=2;DelayUpdate
		// ModifyGraph rgb(wUncageIndicator)=(0,0,0)
		// SetAxis left v_min,(v_max+v_sdev)
		// SetAxis left (:FitResults:y0w[i]-10E-12),(:FitResults:y0w[i]+10E-12)
		SetAxis left (:FitResults:y0w[i]-1.5*abs(wavemin(:FitResults:AmplitudeW))),(:FitResults:y0w[i]+abs(wavemin(:FitResults:AmplitudeW)))

		ModifyGraph width=185,height=80
		ModifyGraph margin=50
		ModifyGraph margin(top)=0,margin(right)=0
		ModifyGraph lblMargin=10

		i += 1
	while (i < vNumStim)






	dowindow/k FitAmplitude
	display /n=FitAmplitude :FitResults:AmplitudeW
	setscale /p x,((numpnts(:FitResults:AmplitudeW)-1)*UncageSpacingV),-UncageSpacingV,"nm",:FitResults:AmplitudeW
	Label bottom "distance \\u#2 (nm)"
	ModifyGraph mode=3,marker=19;DelayUpdate
	//	ErrorBars/T=0/L=0.7 AmplitudeW Y,wave=(AmplitudeSeW,AmplitudeSeW)
	Label left "I\\u#2 (pA)"
	ModifyGraph width=185,height=80
	ModifyGraph margin=50
	ModifyGraph margin(top)=0,margin(right)=0
	ModifyGraph lblMargin=10
	//AutoPositionWindow/M=1/R=graph1

	dowindow/k UncageResponses
	plotrows(:FitResults:ResponseWave2d)
	dowindow /c UncageResponses
	ModifyGraph width=185,height=80
	ModifyGraph margin=50
	ModifyGraph margin(top)=0,margin(right)=0
	ModifyGraph lblMargin=10

	//AutoPositionWindow/M=0/R=graph4 graph5

	dowindow/k UncageFits
	plotrows(:FitResults:ModelPredictionWave2d)
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

	make/o/t/n=(4,2) SummaryFigTableData
	summaryFigTableData[][0] = {"Pockels (mV)","Laser (mW)","ID","notes","OrderIsReversed"}
	summaryFigTableData[][1] = {"","","",""}
	summaryFigTableData[0][1] = {num2str(vPockelsVoltage[0]),num2str(VLaserPower[0]),uid[0],"model fixed kinetics",num2str(UncagingOrderIsReversed)}
	dowindow /k SummaryFigTable
	edit /n=SummaryFigTable SummaryFigTableData
	modifytable autosize = {0,1,-1,0,0}
	ModifyTable width(Point)=0
	modifytable alignment=0

	if(!exists("nResponsePanelColsW"))
		make /n=1 nResponsePanelColsW
		nResponsePanelColsW = 2
	endif
	DoMakeLayout()


	pathinfo $sTemp
	if(v_flag==0)
		if (StringMatch(OutputFormat, "png" ))
			SavePICT/O/E=-5/B=288  /win=SummaryFig as (uid[0] +".png")
		endif
		if (StringMatch(OutputFormat, "pdf"))
			SavePICT/O/E=-2  /win=SummaryFig as (uid[0] +".pdf")
		endif
		if (StringMatch(OutputFormat, "svg" ))
			SavePICT/O/E=-9  /win=SummaryFig as (uid[0] +".svg")
		endif
	else
		if (StringMatch(OutputFormat, "png" ))
			SavePICT/O/E=-5/B=288 /p=$sTemp /win=SummaryFig as (uid[0] +".png")
		endif
		if (StringMatch(OutputFormat, "pdf" ))
			SavePICT/O/E=-2 /p=$sTemp /win=SummaryFig as (uid[0] +".pdf")
		endif
		if (StringMatch(OutputFormat, "svg" ))
			SavePICT/O/E=-9 /p=$sTemp /win=SummaryFig as (uid[0] +".svg")
		endif
		//SavePICT/O/E=-5/B=288 /p=OutputDir /win=SummaryFig as (uid[0] +".png")
	endif


	DoWindow/K SpatialRange
	display /N = SpatialRange :fitresults:amplitudew
	make /o /n=3 w_coef
	K0 = 0;
	CurveFit/H="100"/NTHR=0/TBOX=768 exp_XOffset  :FitResults:AmplitudeW /D
	redimension /n=(1,3) w_coef
	DoWindow /K SummaryInfo
	make /o /n=1 nResponse = numpnts(:FitResults:AmplitudeW)
	edit /N=SummaryInfo uid,vPockelsVoltage,w_coef,nResponse
	AutoPositionWindow /m=0 /r = SummaryFig SummaryInfo
	AutoPositionWindow /m=0 /r = SummaryInfo SpatialRange
	DoWindow /k AmplitudeData
	Edit/K=0 /n=AmplitudeData :FitResults:AmplitudeW
	autopositionwindow /m=1 /r=SummaryInfo AmplitudeData
	
	redimension /n=3 w_coef
	concatenate /o/np {vPockelsVoltage,w_coef,nResponse,:fitresults:amplitudew},w_out
	dowindow /k ExcelOutput
	edit /n=ExcelOutput w_out,uid
endmacro

Macro MakeNormalizedResponse()
	string ResponseDfList = "d0;d100;d200;d300;d400;d500;d600;d700;d800;d900;d1000;d1100;d1200;d1300;d1400"
	string ResponseWaveList = "posNo_0;posNo_1;posNo_2;posNo_3;posNo_4;posNo_5;posNo_6;posNo_7;posNo_8;posNo_9;posNo_10;posNo_11;posNo_12;posNo_13;posNo_14;"
	wavestats :FitResults:AmplitudeW

	duplicate/o :FitResults:y0w y0w

	duplicate /o :UAData:posNo_0 posNo_0
	posNo_0 = (posNo_0 - y0w[0])/abs(V_min)
	duplicate /o :UAData:posNo_1 posNo_1
	posNo_1 = (posNo_1 - y0w[1])/abs(V_min)
	duplicate /o :UAData:posNo_2 posNo_2
	posNo_2 = (posNo_2 - y0w[2])/abs(V_min)
	duplicate /o :UAData:posNo_3 posNo_3
	posNo_3 = (posNo_3 - y0w[3])/abs(V_min)
	duplicate /o :UAData:posNo_4 posNo_4
	posNo_4 = (posNo_4 - y0w[4])/abs(V_min)
	duplicate /o :UAData:posNo_5 posNo_5
	posNo_5 = (posNo_5 - y0w[5])/abs(V_min)
	duplicate /o :UAData:posNo_6 posNo_6
	posNo_6 = (posNo_6 - y0w[6])/abs(V_min)
	duplicate /o :UAData:posNo_7 posNo_7
	posNo_7 = (posNo_7 - y0w[7])/abs(V_min)
	duplicate /o :UAData:posNo_8 posNo_8
	posNo_8 = (posNo_8 - y0w[8])/abs(V_min)
	duplicate /o :UAData:posNo_9 posNo_9
	posNo_9 = (posNo_9 - y0w[9])/abs(V_min)
	duplicate /o :UAData:posNo_10 posNo_10
	posNo_10 = (posNo_10 - y0w[10])/abs(V_min)
	duplicate /o :UAData:posNo_11 posNo_11
	posNo_11 = (posNo_11 - y0w[11])/abs(V_min)
	duplicate /o :UAData:posNo_12 posNo_12
	posNo_12 = (posNo_12 - y0w[12])/abs(V_min)
	duplicate /o :UAData:posNo_13 posNo_13
	posNo_13 = (posNo_13 - y0w[13])/abs(V_min)
	duplicate /o :UAData:posNo_14 posNo_14
	posNo_14 = (posNo_14 - y0w[14])/abs(V_min)
	variable i_ = 0
	string FromThis,ToThis
	do
		FromThis = stringfromlist((V_minrowloc-i_),ResponseWaveList)
		ToThis =	 "root:"+stringfromlist(i_,ResponseDfList)+":"+stringfromlist((V_minrowloc-i_),ResponseWaveList)+"_"+ReplaceString("'", uid[0], "")
		duplicate /o $FromThis $ToThis

		i_ += 1
	while (i_ <= v_minrowloc)	// as long as expression is TRUE
EndMacro

menu "macros"
	"DoUncagingAnalysis/1"
	"DoMakeFigures/2"
	"DoSaveResults/3"
	"Clean_Up/4"
	"SetResponseNan/5"
	"RemoveResponse/6"
	"InsertResponse/7"
	"RemoveSpineData"
	"MakeNormalizedResponse"
end
