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
  string /g OutputPathStr = s_path
	uid = uid + s_filename[0,strlen(s_filename)-5]
	UncagingAnalysis(DataWaveList)

endmacro

//******************************************************************************
//******************************************************************************
//
//******************************************************************************
//******************************************************************************

function UserDefineInitialEstimates(ParametersIn,w_coef,UncageTime,y0timeWindow,Amplitude0window, ResponseMaxTime, DelayToResponseStart,FitStop)
wave ParametersIn, w_coef
variable UncageTime, &y0timeWindow, &Amplitude0window, &ResponseMaxTime, &DelayToResponseStart, &FitStop
//
//graph average response and get user input for initial estimates
//
//graph
dowindow /k review
display /n=review ParametersIn
SetAxis bottom (UncageTime-y0timeWindow),FitStop
SetDrawEnv xcoord= bottom;SetDrawEnv dash= 3;DelayUpdate
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
	String uncaging_response_wave_name	//name of the uncaging response wave
	Prompt uncaging_response_wave_name,"uncaging response wave name",popup,DataWaveList+"Some other wave..."
	String uncaging_power_wave_name	//name of the uncaging power wave
	Prompt uncaging_power_wave_name,"uncaging power wave name",popup,DataWaveList+"Some other wave..."
	Variable n_uncaging_pulses //variable that holds the number of uncaging pulses found in uncaging_power_wave
	Variable fit_range = 0.05	//length of the time window used for fitting double exponential
	Prompt fit_range,"Size of time window for fitting"
	Variable i = 0, j = 0	//dummy variable for iteration control
	Variable decay_time_0 = 0.008	//initial estimate for uncaging response decay time
	Prompt decay_time_0,"response decay time initial estimate"
	Variable rise_time_0 = 0.001	//initial estimate for uncaging response rise time
	Prompt rise_time_0,"response decay time initial estimate"
	Variable DelayToResponseStart = 0
	Variable v_amplitude_0
	Variable Amplitude0window = 0.001
	Prompt Amplitude0window,"Window size for amplitude estimate"
	Variable y0timeWindow = 0.01	//time window before uncaging pulse used to estimate y0
	Prompt y0timeWindow,"Window size for y0 estimate."
	Variable UncageTime	//time of uncaging event for current fit
	// Variable response_max_amplitude	//not used
	Variable ResponseMaxTime	//time of peak amplitude
	Variable user_response	//not used
	// Variable V_FitMaxIters = 100	//not used
	Variable threshold	//threshold value used while examining uncaging_power_wave to determine whether an algorithmically found peak is a true pulse
	// Variable cursor_a, cursor_b	//not used
	Make /O/D w_coef = NaN
	Make /O/D w_sigma = NaN
	Variable V_LevelX = 0	//return value from Igor's built-in findlevel function, used for finding uncaging pulses
	Variable k0, k1, k2, k3, k4	//initial estimates for fit parameters
	Variable fit_start, fit_stop	//variables defining start and stop of the fit window
	Variable V_AbortCode = 0, V_FitError=0 //Igor environment variable indicating fitting abort condition
	Variable v_n_fit_points = Inf// number of points included in fit
	Variable n_false_replicates = 50 //number of "fake responses" to fit, per stimulus
	Variable inter_fit_time //time between end of a fit and start of the next fit for consecutive fits to "real" stimuli
	Variable inter_false_fit_time //time between fake responses
	Variable response_max_time_0, v_delay_to_response_start_0 //intial parameter estimates
	Variable v_amplitude, TotalLengthUncaging
	wave w_Resampled, wColumnMeans
	make /o /n=1 vPockelsVoltage
	variable UserSetPar0
	variable v_flag
	prompt UserSetPar0,"Initial parameter estimates",popup,"Interactive;Default;Auto Guess"
	Variable LaswerPowerWaveScaling = 0.1; Prompt LaswerPowerWaveScaling, "Laser power wave scaling"


	String wave_data_list =  "w_amplitude;w_t0;w_decay_time;w_rise_time;w_uncage_time;w_onset_delay;w_y0;w_amplitude_se;"

// prompt user to select stimulus and response waves from among recently loaded waves
// user is also provided the option "some other wave..."
// if the user indicates that the stimulus/response wave is not lasted, then prompt again, offering a list of all waves in the current data folder
	Variable n_data_waves = itemsinlist(wave_data_list)
	doPrompt "",uncaging_response_wave_name,uncaging_power_wave_name
	if(!(cmpstr(uncaging_response_wave_name, "Some other wave..." )*cmpstr(uncaging_power_wave_name, "Some other wave..." )))
	Prompt uncaging_response_wave_name,"uncaging response wave name",popup,wavelist("*",";","")
	Prompt uncaging_power_wave_name,"uncaging power wave name",popup,wavelist("*",";","")
	doPrompt "",uncaging_response_wave_name,uncaging_power_wave_name
	endif


// set default value for UserSetPar0
	UserSetPar0 = 3

// promp user for starting parameters
	DoPrompt "",fit_range,decay_time_0,rise_time_0,Amplitude0window,y0timeWindow,UserSetPar0,LaswerPowerWaveScaling

// set stimulus and response wave references from chosen names
	wave uncaging_response_wave = $uncaging_response_wave_name
	wave uncaging_power_wave = $uncaging_power_wave_name

// it is unclear why I felt the need to duplicate these waves
duplicate /o uncaging_response_wave w_uncage_response
duplicate /o uncaging_power_wave w_uncage_power

	//before we can make waves to put the fit parameters, we need to know how long to make them
	//loop through uncaging events to count the number of uncaging pulses
	//we find the first two pulses manually before entering the do loop
	//as we treat the first pulse slightly different than the rest, treating the first pulses outside of the loop allows us to avoid using an if-then construct inside the loop
	threshold = 0.8*wavemax(uncaging_power_wave)

//calculate average pockels voltage
duplicate /o uncaging_power_wave temp
temp = threshold < uncaging_power_wave
TotalLengthUncaging = sum(temp)
temp = temp * uncaging_power_wave
vPockelsVoltage[0] = sum(temp)/TotalLengthUncaging

	// find rising edge in power trace greater than defined threshold
	findlevel /Q/EDGE=1 /R= (V_LevelX,) uncaging_power_wave, threshold
	//we find the first pulse manually before entering the do loop, so we start the counter at 1
	Make/O/N=1 w_uncage_time
	w_uncage_time = V_LevelX
	i = 1
	do//do1
		// findlevel may occasionaly find false peaks in the noise during the uncaging pulse, but we can reliably find the time where the voltage drops below the threshold
		findlevel /Q/EDGE=2 /R= (V_LevelX,) uncaging_power_wave, threshold
		//the next findlevel searches again for a rising pulse crossing the threshold
		findlevel /Q/EDGE=1 /R= (V_LevelX,) uncaging_power_wave, threshold
		//if no edge is found, V_LevelX will be set to NaN, this condition indicates that we have passed the last uncaging pulse
		if(numtype(V_LevelX))
			break
		endif
		//at this point in the do loop we have succesfully found an uncaging event, so we increment the counter and save the time
		InsertPoints i, 1, w_uncage_time
		w_uncage_time[i] = V_LevelX
		i = i + 1
		//we never have i>100 in experiments so if i>100 something is wrong with the code and we abort
		if(i>100)//if1
			abort("n>100 uncaging pulses found")
		endif//if1
	while(1)//do1
	n_uncaging_pulses = i
	print "number of uncaging events found: ",n_uncaging_pulses

	//create waves to store parameters
	make /o/n=(n_uncaging_pulses) w_fit_start_pt
	make /o/n=(n_uncaging_pulses) w_fit_start_time
	make /o/n=(n_uncaging_pulses) w_fit_stop_time
	make /o/n=(n_uncaging_pulses) w_amplitude
	make /o/n=(n_uncaging_pulses) w_amplitude_se
	make /o/n=(n_uncaging_pulses) w_t0
	make /o/n=(n_uncaging_pulses) w_decay_time
	make /o/n=(n_uncaging_pulses) w_rise_time
	make /o/n=(n_uncaging_pulses) w_y0
	make /o/n=(n_uncaging_pulses) w_onset_delay
	make /o/n=(n_uncaging_pulses) w_amplitude_1
	make /o/n=(n_uncaging_pulses) w_amplitude_0
	make /o/n=(n_uncaging_pulses) w_amplitude_1_se
	make /o/n=(n_uncaging_pulses,6) w2dFitUncagingResponseCoef
	make /o/n=(n_uncaging_pulses,6) w2dFitUncagingResponseCoefSE
	make /o/n=(n_uncaging_pulses,6) w2dFitUncgRespCoefBootSE

// To get initial estimates of the parameters, the model is first fit to the average response
// to calculate the average response, we need to align the uncaging response relative to the uncaging pulse

//for pre allocating the array to hold the responses, we need to know the number of points, rather than the time
//since the time window is a real number, it's possible that the integral index ranges could have different lengths due to rounding
//if the number of points differs, we use the smallest for analysis
//!note maybe we want to just throw an error if the number of points differs instead
	for(i=0;i < n_uncaging_pulses;i+=1)	// for1
		fit_start = w_uncage_time[i] - y0timeWindow
		w_fit_start_time[i] = fit_start
		UncageTime = fit_start + y0timeWindow
		fit_stop = fit_start + fit_range
		w_fit_stop_time[i] = fit_stop
		v_n_fit_points = min(v_n_fit_points,	(x2pnt(uncaging_response_wave, fit_stop )-x2pnt(uncaging_response_wave, fit_start )))
		w_fit_start_pt[i] = x2pnt(uncaging_response_wave, fit_start )
	endfor		//for1

//copy uncaging responses into 2d wave
make /o/n=(n_uncaging_pulses, v_n_fit_points) w2d_responses
make /o/n=(n_uncaging_pulses, v_n_fit_points) w2d_fits
make /o/n=(v_n_fit_points) w_fit
for(i=0;i < n_uncaging_pulses;i+=1)	// for1
w2d_responses[i][] = uncaging_response_wave[w_fit_start_pt[i]+q]
endfor		//for1
setscale /p y,0,dimdelta(uncaging_response_wave,0),w2d_responses
setscale /p y,0,dimdelta(uncaging_response_wave,0),w2d_fits
setscale /p x,0,dimdelta(uncaging_response_wave,0),w_fit


//calculate average response
make /o /n=(n_uncaging_pulses) w_temp
w_temp = 1
MatrixOp/O w_avg_response=w_temp^t x w2d_responses
redimension /n=(v_n_fit_points) w_avg_response
w_avg_response = w_avg_response/n_uncaging_pulses
setscale /p x,0, dimdelta(uncaging_response_wave,0), w_avg_response

fit_stop = fit_range

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
	UserDefineInitialEstimates(w_avg_response,w_coef,y0timeWindow,y0timeWindow,Amplitude0window, ResponseMaxTime, DelayToResponseStart,fit_stop)
	response_max_time_0 = ResponseMaxTime
	v_delay_to_response_start_0 = DelayToResponseStart
	V_AbortCode = 0
	V_FitError = 0
	FuncFit/N/Q/H="000001" /NTHR=0 DiffTwoExp2 W_coef  w_avg_response /D
	v_delay_to_response_start_0 = w_coef[1]
	v_amplitude_0 = w_coef[0]
	decay_time_0 = w_coef[2]
	rise_time_0 = w_coef[3]
	endif
	if(UserSetPar0 == 2)
	// use parameter values set in dialog
		// ResponseMaxTime = y0timeWindow + 3 * rise_time_0
		DelayToResponseStart = 0
		v_amplitude_0 = -10
    // duplicate /o /r=(0,y0timeWindow) w_avg_response WTemp
    // wavestats /q WTemp
		w_coef = {-10,DelayToResponseStart,rise_time_0,decay_time_0,mean(w_avg_response,0,y0timeWindow),y0timeWindow}
		// perform fit of average response for display purposes
		FuncFit/N/Q/H="000001" /NTHR=0 DiffTwoExp2 W_coef  w_avg_response /D
		ResponseMaxTime = y0timeWindow + v_delay_to_response_start_0 + (decay_time_0*rise_time_0)/(decay_time_0-rise_time_0)*ln(decay_time_0/rise_time_0)
		response_max_time_0 = ResponseMaxTime
		DelayToResponseStart = w_coef[1]
		v_delay_to_response_start_0 = DelayToResponseStart
	endif
	if(UserSetPar0==3)
	// auto guess start parameters
	ResponseMaxTime = y0timeWindow + 3 * rise_time_0
	DelayToResponseStart = 0
	v_amplitude_0 = -10
	w_coef = {-10,DelayToResponseStart,rise_time_0,decay_time_0,mean(w_avg_response,0,y0timeWindow),y0timeWindow}
	FuncFit/N/Q/H="000001" /NTHR=0 DiffTwoExp2 W_coef  w_avg_response /D
	duplicate /o w_coef wAvgUncageResponseFitCoef
	duplicate /o w_sigma wAvgUncageResponseFitCoefSE
	v_delay_to_response_start_0 = w_coef[1]
	DelayToResponseStart = v_delay_to_response_start_0
	v_amplitude_0 = w_coef[0]
	decay_time_0 = w_coef[2]
	rise_time_0 = w_coef[3]
	ResponseMaxTime = y0timeWindow + w_coef[1] + (w_coef[2]*w_coef[3])/(w_coef[2]-w_coef[3])*ln(w_coef[2]/w_coef[3])
	response_max_time_0 = ResponseMaxTime

	endif


	// ============================================================================
	// ============================================================================
	// do fits
	// ============================================================================
	// ============================================================================

Make/O/T/N=2 T_Constraints
T_Constraints[0] = {"K1 > 0","K1 < .01"}

	Prompt user_response, "Is this fit good?", popup, "Yes: Save fit;No: Do a Refit; No response: Save zero; Too noisy, save NaN"

	for(i=0;i < n_uncaging_pulses;i+=1)	// for1
	//the logical control for refitting is handled by performing the fit in a do-while loop with while(user_response=2)
	// set user response = 2 to initially enter the loop
	user_response = 2

// calculate time window for ith uncaging event
		UncageTime = w_uncage_time[i]
		fit_start = UncageTime - y0timeWindow
		fit_stop = fit_start + fit_range

	// set w_coef to initial parameter estimates
		k4 = mean(uncaging_response_wave,fit_start,UncageTime)
		k0 = mean(uncaging_response_wave,(fit_start + response_max_time_0-Amplitude0window),(fit_start + response_max_time_0 + Amplitude0window)) - k4
		k1 = UncageTime + v_delay_to_response_start_0
		k2 = decay_time_0
		k3 = rise_time_0
    W_Coef = {k0, v_delay_to_response_start_0, k2, k3, k4, UncageTime}
		V_AbortCode = 0
		V_FitError = 0

// perform initial fit with some parameters fixed
		// FuncFit/N/Q/H="001101" /NTHR=0 DiffTwoExp2 W_coef  uncaging_response_wave(fit_start, fit_stop) /D /C=T_Constraints
		FuncFit/N/Q/H="011101" /NTHR=0 DiffTwoExp2 W_coef  uncaging_response_wave(fit_start, fit_stop) /D
		// save the amplitude for restricted model
		w_amplitude_1[i] = w_coef[0]
		w_amplitude_1_se[i] = w_sigma[0]
		make /o /n=(v_n_fit_points) w_temp
		w_temp = w2d_responses[i][p]
		setscale /p x,0, dimdelta(uncaging_response_wave,0), w_temp
		duplicate /o w_temp wThisResponse
		//save nonparametric estimate of amplitude
		w_amplitude_0[i] = mean(w_temp,(response_max_time_0-Amplitude0window),(response_max_time_0+Amplitude0window)) - mean(w_temp,0,y0timeWindow)
		// save nonparametric estimate of amplitude
		// ResponseMaxTime = y0timeWindow + w_coef[1] + (w_coef[2]*w_coef[3])/(w_coef[2]-w_coef[3])*ln(w_coef[2]/w_coef[3])
		// w_amplitude_0[i] = mean(uncaging_response_wave,(fit_start+ResponseMaxTime-Amplitude0window),(fit_start+ResponseMaxTime+Amplitude0window)) - mean(uncaging_response_wave,fit_start,(fit_start+y0timeWindow))

do
// open display window for checking the fit; igor will automatically append the fit to the graph when funcfit is called
dowindow /k review
display /n=review uncaging_response_wave[x2pnt(uncaging_response_wave, fit_start ),x2pnt(uncaging_response_wave, fit_stop )]
SetDrawEnv xcoord= bottom;SetDrawEnv dash= 3;DelayUpdate
DrawLine UncageTime,0,UncageTime,1

// perform fit to the full model
// FuncFit/N/Q/H="000001" /NTHR=0 DiffTwoExp2 W_coef  uncaging_response_wave(fit_start, fit_stop) /D /C=T_Constraints
FuncFit/N/Q/H="000001" /NTHR=0 DiffTwoExp2 W_coef  uncaging_response_wave(fit_start, fit_stop) /D


DoUpdate

// allow user to interact to indicate whether the fit is good

user_response = 1 //set default response in the dialog to indicate fit is good
DoPrompt "Goodness of Fit", user_response
if(v_flag)
	abort("user canceled")
endif
switch(user_response)	// numeric switch
case 1:		// fit is good, nothing to do
	break						// exit from switch
case 2: //user indicated to perform a refit, call UserDefineInitialEstimates
	UserDefineInitialEstimates(uncaging_response_wave,w_coef,UncageTime,y0timeWindow,Amplitude0window, ResponseMaxTime, DelayToResponseStart, fit_stop)
	break
case 3: // user indicates no response

	w_coef = 0
	w_sigma = NaN
	break
case 4: //user indicates to save NaN
	w_coef = NaN
	w_sigma = NaN
endswitch

while(user_response == 2) //if user indicated to perform a refit, continue loop, else break

// save fit parameters
		w_amplitude[i] = w_coef[0]
		w_amplitude_se[i] = w_sigma[0]
		w_t0[i] = w_coef[1]
		w_decay_time[i] = w_coef[2]
		w_rise_time[i] = w_coef[3]
		w_y0[i] = w_coef[4]
		w_onset_delay[i] = w_coef[1]
		fit_stop = fit_start + fit_range
		duplicate /o /r=(fit_start, fit_stop) uncaging_response_wave w_t
		w_t = x
		duplicate /o w_t w_fit
		w_fit = DiffTwoExp2(w_coef, w_t)
		w2d_fits[i][] = w_fit[q]
		w2dfituncagingresponsecoef[i][] = w_coef[q]
		w2dFitUncagingResponseCoef[i][] = w_coef[q]
		w2dFitUncagingResponseCoefSE[i][] = w_sigma[q]





	endfor												// for1; loop through uncaging events, performing fits for each

dowindow/k review


return 1
// ============================================================================
// ============================================================================
// perform fits to response trace over periods where no stimulus is given
// ============================================================================
// ============================================================================

make /o /n=((n_uncaging_pulses-1)*n_false_replicates) WNrAmplitude1
make /o /n=((n_uncaging_pulses-1)*n_false_replicates) WNrAmplitude
make /o /n=((n_uncaging_pulses-1)*n_false_replicates) WNrAmplitude0

for(i=0;i < (n_uncaging_pulses-1);i+=1)	// for1
 // there is no stimulus between uncaging pulses
 // loop over uncaging pulses and perform fits during the time between uncaging events
inter_fit_time = w_fit_start_time[i+1] - w_fit_stop_time[i]
inter_false_fit_time = inter_fit_time/n_false_replicates
for(j=0;j < n_false_replicates;j+=1)	// for2
// during each inter-stimulus, perform n_false_replicates fits

// set initial parameters and fit window
fit_start = w_fit_stop_time[i] + j*inter_false_fit_time
UncageTime = fit_start + y0timeWindow
fit_stop = fit_start + fit_range
k4 = mean(uncaging_response_wave,fit_start,UncageTime)
k0 = mean(uncaging_response_wave,(fit_start + response_max_time_0 - Amplitude0window),(fit_start + response_max_time_0 - Amplitude0window)) - k4
v_amplitude = k0
k1 = UncageTime + DelayToResponseStart
k2 = decay_time_0
k3 = rise_time_0
W_Coef = {k0, v_delay_to_response_start_0, k2, k3, k4,UncageTime}
V_AbortCode = 0
V_FitError = 0

// fit restricted model
FuncFit/N/Q/W=2/H="011101" /NTHR=0 DiffTwoExp2 W_coef  uncaging_response_wave(fit_start, fit_stop) /D
// FuncFit/N/Q/W=2/H="001101" /NTHR=0 DiffTwoExp2 W_coef  uncaging_response_wave(fit_start, fit_stop) /D /C=T_Constraints
// ResponseMaxTime = y0timeWindow + w_coef[1] + (w_coef[2]*w_coef[3])/(w_coef[2]-w_coef[3])*ln(w_coef[2]/w_coef[3])
// save nonparametric estimate of amplitude
// WNrAmplitude0[i*n_false_replicates+j] = mean(uncaging_response_wave,(fit_start+ResponseMaxTime-Amplitude0window),(fit_start+ResponseMaxTime+Amplitude0window)) - mean(uncaging_response_wave,fit_start,(fit_start+y0timeWindow))
WNrAmplitude0[i*n_false_replicates+j] = v_amplitude
// save amplitude from restricted model
WNrAmplitude1[i*n_false_replicates+j] = w_coef[0]
// perform full model fit
// FuncFit/N/Q/W=2/H="000001" /NTHR=0 DiffTwoExp2 W_coef  uncaging_response_wave(fit_start, fit_stop) /D /C=T_Constraints
FuncFit/N/Q/W=2/H="000001" /NTHR=0 DiffTwoExp2 W_coef  uncaging_response_wave(fit_start, fit_stop) /D
// duplicate /o/r=[0,4] w_coef w_coef2
// duplicate /o /r =(fit_start, fit_stop) uncaging_response_wave Wtemp
// duplicate /o wTemp wTemp2
// wTemp2 = x
// concatenate /o {wTemp2,wTemp}, ModelFrame
// optimize /q /x=w_coef /R=w_coef /m={0,1} DiffTwoExpSM,w_coef
// w_coef[0,4]=w_coef2[p]
// save amplitude from full model fit
WNrAmplitude[i*n_false_replicates+j] = w_coef[0]
// WNrAmplitude[i*n_false_replicates+j] = w_coef[0]

// save nonparametric amplitude estimate
// WNrAmplitude0[i*n_false_replicates+j] = v_amplitude

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
wave rw2d_response, w_uncage_response,w_uncage_time
wave w_fit_start_time, w_fit_stop_time, w_amplitude, w_amplitude_se, w_t0
wave w_decay_time, w_rise_time, w_y0, w_onset_delay, w_amplitude_1
wave w_amplitude_0, w2d_responses, w2d_fits, rw_uid, WNrAmplitude1
wave WNrAmplitude, WNrAmplitude0, rwPockelsVoltage
wave vPockelsVoltage
variable n_results
string /g OutputPathStr
string /g uid
newpath /o OutputDir, OutputPathStr
SavePICT/O/E=-5/B=288 /p=OutputDir /win=SummaryFig as (uid +".png")
//SavePICT/O/E=-5/B=72 /p=OutputDir /win=SummaryFig as "SummaryFiglr.png"

if(!waveexists(rw2d_response))

make /n=(numpnts(w_uncage_response),8) rw2d_response
make /n=(numpnts(w_uncage_time),8) rw2d_uncage_time
make /n=(numpnts(w_fit_start_time),8) rw2d_fit_start_time
make /n=(numpnts(w_fit_stop_time),8) rw2d_fit_stop_time
make /n=(numpnts(w_amplitude),8) rw2d_fit_amplitude
make /n=(numpnts(w_amplitude_se),8) rw2d_fit_amplitude_se
make /n=(numpnts(w_t0),8) rw2d_fit_t0
make /n=(numpnts(w_decay_time),8) rw2d_fit_decay_time
make /n=(numpnts(w_rise_time),8) rw2d_fit_rise_time
make /n=(numpnts(w_y0),8) rw2d_fit_y0
make /n=(numpnts(w_onset_delay),8) rw2d_fit_onset_delay
make /n=(numpnts(w_amplitude_1),8) rw2d_amplitude_0
make /n=(numpnts(w_amplitude_0),8) rw2d_amplitude_0_np
// make /n=(numpnts(),8)


duplicate w2d_responses rw3d_uncaging_response
redimension /n=(-1,-1,8) rw3d_uncaging_response
duplicate w2d_fits rw3d_fits
redimension /n=(-1,-1,8) rw3d_fits
duplicate WNrAmplitude1 W2dNrAmplitude0
redimension /n=(-1,8) W2dNrAmplitude0
duplicate WNrAmplitude W2dNrAmplitude1
redimension /n=(-1,8) W2dNrAmplitude1
duplicate WNrAmplitude0 W2dNrAmplitude2
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


rw2d_response[][n_results] = w_uncage_response[p]
rw2d_uncage_time[][n_results] = w_uncage_time[p]
rw2d_fit_start_time[][n_results] = w_fit_start_time[p]
rw2d_fit_stop_time[][n_results] = w_fit_stop_time[p]
rw2d_fit_amplitude[][n_results] = w_amplitude[p]
rw2d_fit_amplitude_se[][n_results] = w_amplitude_se[p]
rw2d_fit_t0[][n_results] = w_t0[p]
rw2d_fit_decay_time[][n_results] = w_decay_time[p]
rw2d_fit_rise_time[][n_results] = w_rise_time[p]
rw2d_fit_y0[][n_results] = w_y0[p]
rw2d_fit_onset_delay[][n_results] = w_onset_delay[p]
rw2d_amplitude_0[][n_results] = w_amplitude_1[p]
rw2d_amplitude_0_np[][n_results] = w_amplitude_0[p]
W2dNrAmplitude0[][n_results] = WNrAmplitude1[p]
W2dNrAmplitude1[][n_results] = WNrAmplitude[p]
W2dNrAmplitude2[][n_results] = WNrAmplitude0[p]
// [][n_results] = [p]

rw3d_uncaging_response[][][n_results]=w2d_responses[p][q]
rw3d_fits[][][n_results]=w2d_fits[p][q]

InsertPoints numpnts(WAmplitudeCorrelation), 1, WAmplitudeCorrelation
InsertPoints numpnts(WAmplitudeNrCorrelation), 1, WAmplitudeNrCorrelation
InsertPoints numpnts(rwPockelsVoltage), 1, rwPockelsVoltage
rwPockelsVoltage[n_results] = vPockelsVoltage[0]
WAmplitudeCorrelation[n_results] = StatsCorrelation(w_amplitude,w_amplitude_0)
WAmplitudeNrCorrelation[n_results] = StatsCorrelation(WNrAmplitude,WNrAmplitude0)
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
kill_wave_list("w_uncage_response;w_uncage_power;")
kill_wave_list("fit_w2d_responses;w_t;w2d_fake_pars;fit_w_response_out;w_bs_amp0;w_bs_amp0alt;w_bs_amp0_Hist;w_bs_amp0alt_Hist;")
kill_wave_list("w2d_responses;w2d_fits;w_fit;w_temp;w_avg_response;T_Constraints;fit_w_avg_response;")
kill_wave_list("w_rise_time;w_y0;w_onset_delay;w_amplitude_1;w_amplitude_0;w_amplitude_1_se;")
kill_wave_list("w_fit_start_time;w_fit_stop_time;w_amplitude;w_amplitude_se;w_t0;w_decay_time;")
kill_wave_list("w_response_out;w_power_out;w_coef;W_sigma;w_uncage_time;w_fit_start_pt;")
kill_wave_list("ACH_1;ACH_3;w_uncage_time;w_refs;w_stim1;w_stim2;w_stim3;w_stim4;w_stim5;w_response_out;w_power_out;w2d_responses;w2d_stim;w_temp;w_avg_response;w_avg_power;")
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
wave w_amplitude,w_amplitude_se,w_t0,w_decay_time,w_rise_time,w_y0,w_onset_delay
wave w_amplitude_1,w_amplitude_0,w_amplitude_1_se,w2d_responses,w2d_responses,w2d_fits
variable i_
prompt i_,"Point number"
doprompt "Enter value",i_
//ShowInfo/CP=0
//cursor a,$StringFromList(0, tracenamelist("",";",1) ),0
w_amplitude[i_]=NaN;w_amplitude_se[i_]=NaN;w_t0[i_]=NaN;w_decay_time[i_]=NaN;w_rise_time[i_]=NaN;w_y0[i_]=NaN;
w_onset_delay[i_]=NaN;w_amplitude_1[i_]=NaN;w_amplitude_0[i_]=NaN;w_amplitude_1_se[i_]=NaN;
w2d_responses[i_][]=NaN
w2d_fits[i_][]=NaN
end

//******************************************************************************
//******************************************************************************
//
//******************************************************************************
//******************************************************************************

function RemoveResponse()
wave w_amplitude,w_amplitude_se,w_t0,w_decay_time,w_rise_time,w_y0,w_onset_delay
wave w_amplitude_1,w_amplitude_0,w_amplitude_1_se,w2d_responses,w2d_responses,w2d_fits
variable i_
prompt i_,"Point number"
doprompt "Enter value",i_
DeletePoints i_,1, w_amplitude,w_amplitude_se,w_t0,w_decay_time,w_rise_time,w_y0;DelayUpdate
DeletePoints i_,1, w_onset_delay,w_amplitude_1,w_amplitude_0,w_amplitude_1_se
DeletePoints i_,1, w2d_responses
DeletePoints i_,1, w2d_fits
end

//******************************************************************************
//******************************************************************************
//
//******************************************************************************
//******************************************************************************

function InsertResponse()
wave w_amplitude,w_amplitude_se,w_t0,w_decay_time,w_rise_time,w_y0,w_onset_delay
wave w_amplitude_1,w_amplitude_0,w_amplitude_1_se,w2d_responses,w2d_responses,w2d_fits
Insertpoints 0,1, w_amplitude,w_amplitude_se,w_t0,w_decay_time,w_rise_time,w_y0;DelayUpdate
Insertpoints 0,1, w_onset_delay,w_amplitude_1,w_amplitude_0,w_amplitude_1_se
Insertpoints 0,1, w2d_responses
Insertpoints 0,1, w2d_fits

w_amplitude[0]=NaN;w_amplitude_se[0]=NaN;w_t0[0]=NaN;w_decay_time[0]=NaN;w_rise_time[0]=NaN;w_y0[0]=NaN;
w_onset_delay[0]=NaN;w_amplitude_1[0]=NaN;w_amplitude_0[0]=NaN;w_amplitude_1_se[0]=NaN;
w2d_responses[0][]=NaN
w2d_fits[0][]=NaN

end

//******************************************************************************
//******************************************************************************
//
//******************************************************************************
//******************************************************************************

macro DoMakeFigures()
	//MakeFigures()

	variable vNumStim
	variable vUncageSpacing = 400
	variable i
	String cmd

	//vPockelsVoltage = 20*vPockelsVoltage
	dowindow/k FullResponse
	display /n=FullResponse w_uncage_response
	Label bottom "time \\u#2 (s)"
	Label left "I\\u#2 (pA)"
	wavestats /	q w_uncage_response
	SetAxis left v_min,(v_max+v_sdev)
	duplicate /o w_uncage_time wUncageIndicator
	wUncageIndicator = (v_max+v_sdev)
	appendtograph wUncageIndicator vs w_uncage_time
	ModifyGraph mode(wUncageIndicator)=3,marker(wUncageIndicator)=2;DelayUpdate
	ModifyGraph rgb(wUncageIndicator)=(0,0,0)
	ModifyGraph width={Aspect,7}
	ModifyGraph width=430,height=225
	ModifyGraph margin(bottom)=40
	ModifyGraph margin(left)=40
	ModifyGraph margin(top)=20
	ModifyGraph margin(right)=20

	vNumStim = dimsize(w2d_responses,0)
	duplicate /o w_uncage_time wUncageIndicatorTime
	wUncageIndicatorTime = w_uncage_time[0]-w_fit_start_time[0]
	i = 0
	do
	sprintf cmd, "dowindow /k response%s", num2str(i)
	print cmd
	Execute cmd
	sprintf cmd, "Display /n= response%s w2d_responses[%s][*]", num2str(i), num2str(i)
	print cmd
	Execute cmd
		Label bottom "\\u#2time (ms)"
		Label left "I\\u#2 (pA)"
		appendtograph w2d_fits[i][*]
		ModifyGraph rgb(w2d_fits)=(0,0,0)
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
	display /n=FitAmplitude w_amplitude
	setscale /p x,((numpnts(w_amplitude)-1)*vUncageSpacing),-vUncageSpacing,"nm",w_amplitude
	Label bottom "distance \\u#2 (nm)"
	ModifyGraph mode=3,marker=19;DelayUpdate
//	ErrorBars/T=0/L=0.7 w_amplitude Y,wave=(w_amplitude_se,w_amplitude_se)
	Label left "I\\u#2 (pA)"
	SetAxis left *,max(wavemax(w_amplitude),0)
	Sort WNrAmplitude WNrAmplitude
	setscale /i x,0,1,WNrAmplitude
	SetDrawEnv ycoord= left;SetDrawEnv dash= 3;DelayUpdate
	DrawLine 0,(WNrAmplitude(0.01)),1,(WNrAmplitude(0.01))
	SetDrawEnv ycoord= left;SetDrawEnv dash= 3;DelayUpdate
	DrawLine 0,(WNrAmplitude(0.05)),1,(WNrAmplitude(0.05))
	SetDrawEnv ycoord= left;SetDrawEnv dash= 3;DelayUpdate
	DrawLine 0,(WNrAmplitude(0.5)),1,(WNrAmplitude(0.5))
	ModifyGraph width=185,height=80
	ModifyGraph margin=50
	ModifyGraph margin(top)=0,margin(right)=0
	ModifyGraph lblMargin=10
	//AutoPositionWindow/M=1/R=graph1

	dowindow/k UncageResponses
	plotrows(w2d_responses)
	dowindow /c UncageResponses
	ModifyGraph width=185,height=80
	ModifyGraph margin=50
	ModifyGraph margin(top)=0,margin(right)=0
	ModifyGraph lblMargin=10

	//AutoPositionWindow/M=0/R=graph4 graph5

	dowindow/k UncageFits
	plotrows(w2d_fits)
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



endmacro


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
macro MakeLayout()
variable i,nResponsePanelCols,j
string cmd

nResponsePanelCols = 3

dowindow /k SummaryFig
newlayout /n=SummaryFig
ModifyGraph /w=FullResponse width=(245*nResponsePanelCols-60)
appendtolayout FullResponse
modifylayout left(FullResponse)=0
modifylayout top(FullResponse)=0

duplicate /o w_amplitude UncagePosition
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


if(i*nResponsePanelCols + j + 1 >= dimsize(w2d_fits,0))
break
endif

i += 1
while(i < dimsize(w2d_fits,0))

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
