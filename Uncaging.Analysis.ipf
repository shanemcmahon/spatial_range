#pragma rtGlobals=3		// Use modern global access method and strict wave access.

Function User_Continue(ctrlName) : ButtonControl
	String ctrlName
	DoWindow/K pause_for_user
End

macro Do_Uncaging_Analysis()
// helper macro for uncaging_analysis Function
// the global variables s_wavenames and s_path are not available inside functions, therefore somethings are done first in the macro
	String /g data_wave_list = s_wavenames
  String /g uid =  StringFromList(ItemsInList(s_path, ":")-4, s_path , ":") + "_" + StringFromList(ItemsInList(s_path, ":")-3, s_path , ":") + "_" + StringFromList(ItemsInList(s_path, ":")-2, s_path , ":")+"_ts"+s_path[strlen(s_path)-4,strlen(s_path)-2]
	Uncaging_Analysis(data_wave_list)

endmacro

function User_Define_Initial_Estimates(w_in,w_coef,uncage_time,y0_time_window,v_amplitude_0_window, response_max_time, v_delay_to_response_start,FitStop)
wave w_in, w_coef
variable uncage_time, &y0_time_window, &v_amplitude_0_window, &response_max_time, &v_delay_to_response_start, &FitStop
//
//graph average response and get user input for initial estimates
//
//graph
dowindow /k review
display /n=review w_in
SetAxis bottom (uncage_time-y0_time_window),FitStop
SetDrawEnv xcoord= bottom;SetDrawEnv dash= 3;DelayUpdate
//draw line indicates the location of the uncaging pulse in the aligned average
//because we have aligned the resonse traces at the begining of the fit window, the uncaging pulse occurs at time = uncage_time
DrawLine uncage_time,0,uncage_time,1
ShowInfo/CP=0/W=review
cursor a,$StringFromList(0, tracenamelist("",";",1) ),(uncage_time-y0_time_window)
//user interaction
//estimate response onset delay
NewPanel/K=2 /n=pause_for_user as "Pause for user"; AutoPositionWindow/M=1/R=review
DrawText 21,20,"Adjust cursor A to estimate response start";	DrawText 21,40,"time."
Button button0,pos={80,58},size={92,20},title="Continue"; Button button0,proc=User_Continue
PauseForUser pause_for_user, review
v_delay_to_response_start = xcsr(a)-uncage_time
k1 = xcsr(a)
//estimate peak location and amplitude
NewPanel/K=2 /n=pause_for_user as "Pause for user"; AutoPositionWindow/M=1/R=review
DrawText 21,20,"Adjust cursor A to peak response.";//	DrawText 21,40,"Click Continue."
Button button0,pos={80,58},size={92,20},title="Continue"; Button button0,proc=User_Continue
PauseForUser pause_for_user, review
response_max_time = xcsr(a)
k0 = mean(w_in,response_max_time-v_amplitude_0_window,response_max_time+v_amplitude_0_window) - mean(w_in,0,uncage_time)
k3 = (response_max_time-k1)*0.33
//estimate decay time
NewPanel/K=2 /n=pause_for_user as "Pause for user"; AutoPositionWindow/M=1/R=review
DrawText 21,20,"Adjust cursor A to indicate time at which the";	DrawText 21,40," response has decayed by 90%."
Button button0,pos={80,58},size={92,20},title="Continue"; Button button0,proc=User_Continue
PauseForUser pause_for_user, review
k2 = (xcsr(a)-response_max_time)*0.33

k4 = mean(w_in,0,uncage_time)
W_Coef = {k0, v_delay_to_response_start, k2, k3, k4, uncage_time}

cursor a,$StringFromList(0, tracenamelist("",";",1) ),(FitStop)


NewPanel/K=2 /n=pause_for_user as "Pause for user"; AutoPositionWindow/M=1/R=review
DrawText 21,20,"If necessary, adjust cursor A to edit fit range";	DrawText 21,40,""
Button button0,pos={80,58},size={92,20},title="Continue"; Button button0,proc=User_Continue
PauseForUser pause_for_user, review
FitStop = xcsr(a)

end

function Uncaging_Analysis(data_wave_list)
	String data_wave_list
	String uncaging_response_wave_name	//name of the uncaging response wave
	Prompt uncaging_response_wave_name,"uncaging response wave name",popup,data_wave_list+"Some other wave..."
	String uncaging_power_wave_name	//name of the uncaging power wave
	Prompt uncaging_power_wave_name,"uncaging power wave name",popup,data_wave_list+"Some other wave..."
	Variable n_uncaging_pulses //variable that holds the number of uncaging pulses found in uncaging_power_wave
	Variable fit_range = 0.05	//length of the time window used for fitting double exponential
	Prompt fit_range,"Size of time window for fitting"
	Variable i = 0, j = 0	//dummy variable for iteration control
	Variable decay_time_0 = 0.008	//initial estimate for uncaging response decay time
	Prompt decay_time_0,"response decay time initial estimate"
	Variable rise_time_0 = 0.001	//initial estimate for uncaging response rise time
	Prompt rise_time_0,"response decay time initial estimate"
	Variable v_delay_to_response_start = 0
	Variable v_amplitude_0
	Variable v_amplitude_0_window = 0.001
	Prompt v_amplitude_0_window,"Window size for amplitude estimate"
	Variable y0_time_window = 0.01	//time window before uncaging pulse used to estimate y0
	Prompt y0_time_window,"Window size for y0 estimate."
	Variable uncage_time	//time of uncaging event for current fit
	// Variable response_max_amplitude	//not used
	Variable response_max_time	//time of peak amplitude
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
	prompt UserSetPar0,"Initial parameter estimates",popup,"Interactive;Default;Auto Guess"


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
	DoPrompt "",fit_range,decay_time_0,rise_time_0,v_amplitude_0_window,y0_time_window,UserSetPar0

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
		fit_start = w_uncage_time[i] - y0_time_window
		w_fit_start_time[i] = fit_start
		uncage_time = fit_start + y0_time_window
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
	// to get user input for initial paramters we call User_Define_Initial_Estimates
	// the function has inputs (w_in,w_coef,uncage_time,y0_time_window,v_amplitude_0_window)
	// //because we have aligned the resonse traces at the begining of the fit window, the uncaging pulse occurs at time = y0_time_window
	User_Define_Initial_Estimates(w_avg_response,w_coef,y0_time_window,y0_time_window,v_amplitude_0_window, response_max_time, v_delay_to_response_start,fit_stop)
	response_max_time_0 = response_max_time
	v_delay_to_response_start_0 = v_delay_to_response_start
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
		// response_max_time = y0_time_window + 3 * rise_time_0
		v_delay_to_response_start = 0
		v_amplitude_0 = -10
    // duplicate /o /r=(0,y0_time_window) w_avg_response WTemp
    // wavestats /q WTemp
		w_coef = {-10,v_delay_to_response_start,rise_time_0,decay_time_0,mean(w_avg_response,0,y0_time_window),y0_time_window}
		// perform fit of average response for display purposes
		FuncFit/N/Q/H="000001" /NTHR=0 DiffTwoExp2 W_coef  w_avg_response /D
		response_max_time = y0_time_window + v_delay_to_response_start_0 + (decay_time_0*rise_time_0)/(decay_time_0-rise_time_0)*ln(decay_time_0/rise_time_0)
		response_max_time_0 = response_max_time
		v_delay_to_response_start = w_coef[1]
		v_delay_to_response_start_0 = v_delay_to_response_start
	endif
	if(UserSetPar0==3)
	// auto guess start parameters
	response_max_time = y0_time_window + 3 * rise_time_0
	v_delay_to_response_start = 0
	v_amplitude_0 = -10
	w_coef = {-10,v_delay_to_response_start,rise_time_0,decay_time_0,mean(w_avg_response,0,y0_time_window),y0_time_window}
	FuncFit/N/Q/H="000001" /NTHR=0 DiffTwoExp2 W_coef  w_avg_response /D
	duplicate /o w_coef wAvgUncageResponseFitCoef
	duplicate /o w_sigma wAvgUncageResponseFitCoefSE
	v_delay_to_response_start_0 = w_coef[1]
	v_delay_to_response_start = v_delay_to_response_start_0
	v_amplitude_0 = w_coef[0]
	decay_time_0 = w_coef[2]
	rise_time_0 = w_coef[3]
	response_max_time = y0_time_window + w_coef[1] + (w_coef[2]*w_coef[3])/(w_coef[2]-w_coef[3])*ln(w_coef[2]/w_coef[3])
	response_max_time_0 = response_max_time

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
		uncage_time = w_uncage_time[i]
		fit_start = uncage_time - y0_time_window
		fit_stop = fit_start + fit_range

	// set w_coef to initial parameter estimates
		k4 = mean(uncaging_response_wave,fit_start,uncage_time)
		k0 = mean(uncaging_response_wave,(fit_start + response_max_time_0-v_amplitude_0_window),(fit_start + response_max_time_0 + v_amplitude_0_window)) - k4
		k1 = uncage_time + v_delay_to_response_start_0
		k2 = decay_time_0
		k3 = rise_time_0
    W_Coef = {k0, v_delay_to_response_start_0, k2, k3, k4, uncage_time}
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
		w_amplitude_0[i] = mean(w_temp,(response_max_time_0-v_amplitude_0_window),(response_max_time_0+v_amplitude_0_window)) - mean(w_temp,0,y0_time_window)
		// save nonparametric estimate of amplitude
		// response_max_time = y0_time_window + w_coef[1] + (w_coef[2]*w_coef[3])/(w_coef[2]-w_coef[3])*ln(w_coef[2]/w_coef[3])
		// w_amplitude_0[i] = mean(uncaging_response_wave,(fit_start+response_max_time-v_amplitude_0_window),(fit_start+response_max_time+v_amplitude_0_window)) - mean(uncaging_response_wave,fit_start,(fit_start+y0_time_window))

do
// open display window for checking the fit; igor will automatically append the fit to the graph when funcfit is called
dowindow /k review
display /n=review uncaging_response_wave[x2pnt(uncaging_response_wave, fit_start ),x2pnt(uncaging_response_wave, fit_stop )]
SetDrawEnv xcoord= bottom;SetDrawEnv dash= 3;DelayUpdate
DrawLine uncage_time,0,uncage_time,1

// perform fit to the full model
// FuncFit/N/Q/H="000001" /NTHR=0 DiffTwoExp2 W_coef  uncaging_response_wave(fit_start, fit_stop) /D /C=T_Constraints
FuncFit/N/Q/H="000001" /NTHR=0 DiffTwoExp2 W_coef  uncaging_response_wave(fit_start, fit_stop) /D


DoUpdate

// allow user to interact to indicate whether the fit is good

user_response = 1 //set default response in the dialog to indicate fit is good
DoPrompt "Goodness of Fit", user_response
switch(user_response)	// numeric switch
case 1:		// fit is good, nothing to do
	break						// exit from switch
case 2: //user indicated to perform a refit, call User_Define_Initial_Estimates
	User_Define_Initial_Estimates(uncaging_response_wave,w_coef,uncage_time,y0_time_window,v_amplitude_0_window, response_max_time, v_delay_to_response_start, fit_stop)
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


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
// experimental code for bootstrap standard errors
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
// 		duplicate /o w_coef w_coef0
// 		duplicate /o wThisResponse wThisResponseResiduals
// 		wThisResponseResiduals = wThisResponse - w_fit
// 		duplicate /o wThisResponse wThisResponse2
// 		w_coef0[5] = y0_time_window
// 		variable vNBootStrapResamples = 100
// 		make /o /n=(vNBootStrapResamples,6) wCoefBoot
// 	for(j=0;j<vNBootStrapResamples;j += 1)	//for2, bootstrap
// 	statsresample /N =(v_n_fit_points) wThisResponseResiduals
// 	wThisResponse2 = wThisResponse + w_Resampled
// 	w_coef = w_coef0
// 	FuncFit/N/Q/w=2/H="000001" /NTHR=0 DiffTwoExp2 W_coef  wThisResponse2 /D
// 	wCoefBoot[j][] = w_coef[q]
// 	endfor		//for2
// ColumnMeans(wCoefBoot)
// duplicate /o wCoefBoot wCoefBootDiffs
// wCoefBootDiffs = wCoefBoot - wColumnMeans[q]
// wCoefBootDiffs = wCoefBootDiffs^2
// ColumnMeans(wCoefBootDiffs)
// wColumnMeans = wColumnMeans^0.5
// w2dFitUncgRespCoefBootSE[i][] = wColumnMeans[q]


	endfor												// for1; loop through uncaging events, performing fits for each

dowindow/k review



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
uncage_time = fit_start + y0_time_window
fit_stop = fit_start + fit_range
k4 = mean(uncaging_response_wave,fit_start,uncage_time)
k0 = mean(uncaging_response_wave,(fit_start + response_max_time_0 - v_amplitude_0_window),(fit_start + response_max_time_0 - v_amplitude_0_window)) - k4
v_amplitude = k0
k1 = uncage_time + v_delay_to_response_start
k2 = decay_time_0
k3 = rise_time_0
W_Coef = {k0, v_delay_to_response_start_0, k2, k3, k4,uncage_time}
V_AbortCode = 0
V_FitError = 0

// fit restricted model
FuncFit/N/Q/W=2/H="011101" /NTHR=0 DiffTwoExp2 W_coef  uncaging_response_wave(fit_start, fit_stop) /D
// FuncFit/N/Q/W=2/H="001101" /NTHR=0 DiffTwoExp2 W_coef  uncaging_response_wave(fit_start, fit_stop) /D /C=T_Constraints
// response_max_time = y0_time_window + w_coef[1] + (w_coef[2]*w_coef[3])/(w_coef[2]-w_coef[3])*ln(w_coef[2]/w_coef[3])
// save nonparametric estimate of amplitude
// WNrAmplitude0[i*n_false_replicates+j] = mean(uncaging_response_wave,(fit_start+response_max_time-v_amplitude_0_window),(fit_start+response_max_time+v_amplitude_0_window)) - mean(uncaging_response_wave,fit_start,(fit_start+y0_time_window))
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

macro DoMakeFigures()
MakeFigures()

dowindow/k graph0
display w_amplitude
setscale /p x,((numpnts(w_amplitude)-1)*150),-150,"nm",w_amplitude
Label bottom "distance \\u#2 (nm)"
ModifyGraph mode=3,marker=19;DelayUpdate
ErrorBars w_amplitude Y,wave=(w_amplitude_se,w_amplitude_se)
SetAxis left wavemin(w_amplitude),wavemax(w_amplitude)
Label left "I\\u#2 (pA)"
Sort WNrAmplitude WNrAmplitude
setscale /i x,0,1,WNrAmplitude
SetDrawEnv ycoord= left;SetDrawEnv dash= 3;DelayUpdate
DrawLine 0,(WNrAmplitude(0.01)),1,(WNrAmplitude(0.01))
SetAxis left *,0



endmacro

Function MakeFigures()
wave w2dFitUncagingResponseCoef, w_amplitude_0,w2dFitUncagingResponseCoefSE,w2d_fits
wave w2d_responses,w_amplitude,w_amplitude_1,w_amplitude_1_se,w_amplitude_se
wave w_decay_time,w_fit_start_pt,w_fit_start_time,w_fit_stop_time,w_onset_delay
wave w_rise_time, w_t0,w_uncage_time,w_y0
variable VXPos,UserResponse
variable VNCols = DimSize(w2dFitUncagingResponseCoef, 0 )
Prompt UserResponse, "Do you whish to set any points to NaN?", popup, "Select point for deletion; Continue analysis"
dowindow/k graph0
display w_amplitude
ModifyGraph mode=3,marker=19;
SetDrawEnv ycoord= left;SetDrawEnv dash= 3;
showinfo
cursor a,$StringFromList(0, tracenamelist("",";",1) ),0
UserResponse = 1

do
DoUpdate
doprompt "",UserResponse
if(UserResponse!=1)
break
endif

NewPanel/K=2 /n=pause_for_user as "Pause for user"; AutoPositionWindow/M=0/R=graph0
DrawText 21,20,"Select points for deletion with cursor a";	DrawText 21,40,""
Button button0,pos={80,58},size={92,20},title="Continue"; Button button0,proc=User_Continue
PauseForUser pause_for_user, graph0

VXPos = xcsr(a)
if(VXPos==(VNCols-1))
DeletePoints VXPos,1, w2dFitUncagingResponseCoef
InsertPoints 0,1, w2dFitUncagingResponseCoef
DeletePoints VXPos,1,w2dFitUncagingResponseCoefSE
InsertPoints 0,1,w2dFitUncagingResponseCoefSE
DeletePoints VXPos,1,w2d_fits
InsertPoints 0,1,w2d_fits
DeletePoints VXPos,1,w2d_responses
InsertPoints 0,1,w2d_responses

w2dFitUncagingResponseCoef[0][] = NaN
w2dFitUncagingResponseCoefSE[0][] = NaN
w2d_fits[0][] = NaN
w2d_responses[0][] = NaN

DeletePoints VXPos,1,w_amplitude
InsertPoints 0,1,w_amplitude
DeletePoints VXPos,1, w_amplitude_0
InsertPoints 0,1, w_amplitude_0
DeletePoints VXPos,1, w_amplitude_1
InsertPoints 0,1, w_amplitude_1
DeletePoints VXPos,1, w_amplitude_1_se
InsertPoints 0,1, w_amplitude_1_se
DeletePoints VXPos,1, w_amplitude_se
InsertPoints 0,1, w_amplitude_se
DeletePoints VXPos,1, w_decay_time
InsertPoints 0,1, w_decay_time
DeletePoints VXPos,1, w_fit_start_pt
InsertPoints 0,1, w_fit_start_pt
DeletePoints VXPos,1, w_fit_start_time
InsertPoints 0,1, w_fit_start_time
DeletePoints VXPos,1, w_fit_stop_time
InsertPoints 0,1, w_fit_stop_time
DeletePoints VXPos,1, w_onset_delay
InsertPoints 0,1, w_onset_delay
DeletePoints VXPos,1, w_rise_time
InsertPoints 0,1, w_rise_time
DeletePoints VXPos,1, w_t0
InsertPoints 0,1, w_t0
DeletePoints VXPos,1, w_uncage_time
InsertPoints 0,1, w_uncage_time
DeletePoints VXPos,1, w_y0
InsertPoints 0,1, w_y0

w_amplitude[0] = NaN
w_amplitude_0[0] = NaN
w_amplitude_1[0] = NaN
w_amplitude_1_se[0] = NaN
w_amplitude_se[0] = NaN
w_decay_time[0] = NaN
w_fit_start_pt[0] = NaN
w_fit_start_time[0] = NaN
w_fit_stop_time[0] = NaN
w_onset_delay[0] = NaN
w_rise_time[0] = NaN
w_t0[0] = NaN
w_uncage_time[0] = NaN
w_y0[0] = NaN
else

w2dFitUncagingResponseCoef[(VXPos)][] = NaN
w2dFitUncagingResponseCoefSE[(VXPos)][] = NaN
w2d_fits[(VXPos)][] = NaN
w2d_responses[(VXPos)][] = NaN

w_amplitude[(VXPos)] = NaN
w_amplitude_0[(VXPos)] = NaN
w_amplitude_1[(VXPos)] = NaN
w_amplitude_1_se[(VXPos)] = NaN
w_amplitude_se[(VXPos)] = NaN
w_decay_time[(VXPos)] = NaN
w_fit_start_pt[(VXPos)] = NaN
w_fit_start_time[(VXPos)] = NaN
w_fit_stop_time[(VXPos)] = NaN
w_onset_delay[(VXPos)] = NaN
w_rise_time[(VXPos)] = NaN
w_t0[(VXPos)] = NaN
w_uncage_time[(VXPos)] = NaN
w_y0[(VXPos)] = NaN

endif


while(UserResponse == 1)

end

macro do_save_results()
if(!exists("rw_uid"))
make /o /n=0 /t rw_uid
endif
InsertPoints numpnts(rw_uid), 1, rw_uid
rw_uid[numpnts(rw_uid)] = uid
save_results()
endmacro

function save_results()
wave rw2d_response, w_uncage_response,w_uncage_time
wave w_fit_start_time, w_fit_stop_time, w_amplitude, w_amplitude_se, w_t0
wave w_decay_time, w_rise_time, w_y0, w_onset_delay, w_amplitude_1
wave w_amplitude_0, w2d_responses, w2d_fits, rw_uid, WNrAmplitude1
wave WNrAmplitude, WNrAmplitude0, rwPockelsVoltage
wave vPockelsVoltage
variable n_results

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

macro Kill_Input_Waves()
	//
	//Kill_Input_Waves kills the waves listed in the data_wave_list String
	//if data_wave_list does not exist, then the function will attempt to kill the waves listed in s_wavenames
	//s_wavenames is auto generated by igor whenever a data wave is loaded
	//this function would normally be invoked by the user after Do_Uncaging_Analysis which would kill only the input waves but leave the output from the analysis in place
	//

	if(!exists("data_wave_list"))//if1
		if(!exists("s_wavenames"))//if2
			return
		endif//endif2
		String data_wave_list = s_wavenames
	endif//endif1
	Kill_Wave_List(data_wave_list)
endmacro








menu "macros"
	"Do_Uncaging_Analysis/1"
	"DoMakeFigures/2"
	"Do_Save_Results/3"
	"Clean_Up/4"
	"Kill_Input_Waves"
end
