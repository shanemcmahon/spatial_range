#pragma rtGlobals=3		// Use modern global access method and strict wave access.

Function User_Continue(ctrlName) : ButtonControl
	String ctrlName
	DoWindow/K pause_for_user
End
macro Do_Uncaging_Analysis()
	String /g data_wave_list = s_wavenames
  String /g uid =  StringFromList(ItemsInList(s_path, ":")-4, s_path , ":") + "_" + StringFromList(ItemsInList(s_path, ":")-3, s_path , ":") + "_" + StringFromList(ItemsInList(s_path, ":")-2, s_path , ":")
	Uncaging_Analysis(data_wave_list)

	//Kill_Wave_List(data_wave_list)
endmacro

function User_Define_Initial_Estimates(w_in,w_coef,uncage_time,y0_time_window,v_amplitude_0_window, response_max_time, v_delay_to_response_start)
wave w_in, w_coef
variable uncage_time, &y0_time_window, &v_amplitude_0_window, &response_max_time, &v_delay_to_response_start
//
//graph average response and get user input for initial estimates
//
//graph
dowindow /k review
display /n=review w_in
SetDrawEnv xcoord= bottom;SetDrawEnv dash= 3;DelayUpdate
//draw line indicates the location of the uncaging pulse in the aligned average
//because we have aligned the resonse traces at the begining of the fit window, the uncaging pulse occurs at time = uncage_time
DrawLine uncage_time,0,uncage_time,1
//cursor a,$"w_in",0
ShowInfo/CP=0/W=review

//user interaction
//estimate response onset delay
NewPanel/K=2 /n=pause_for_user as "Pause for user"; AutoPositionWindow/M=1/R=review
DrawText 21,20,"Adjust cursor A to estimate response start";	DrawText 21,40,"time."
Button button0,pos={80,58},size={92,20},title="Continue"; Button button0,proc=User_Continue
PauseForUser pause_for_user, review
v_delay_to_response_start = xcsr(a)-uncage_time
k1 = xcsr(a)-uncage_time

//estimate peak location and amplitude
NewPanel/K=2 /n=pause_for_user as "Pause for user"; AutoPositionWindow/M=1/R=review
DrawText 21,20,"Adjust cursor A to peak response.";//	DrawText 21,40,"Click Continue."
Button button0,pos={80,58},size={92,20},title="Continue"; Button button0,proc=User_Continue
PauseForUser pause_for_user, review
response_max_time = xcsr(a)
// v_amplitude_0 = mean(w_in,response_max_time-v_amplitude_0_window,response_max_time+v_amplitude_0_window) - mean(w_in,0,uncage_time)
k0 = mean(w_in,response_max_time-v_amplitude_0_window,response_max_time+v_amplitude_0_window) - mean(w_in,0,uncage_time)
// rise_time_0 = (response_max_time-v_delay_to_response_start)*0.33
k3 = (response_max_time-v_delay_to_response_start)*0.33
//estimate decay time
NewPanel/K=2 /n=pause_for_user as "Pause for user"; AutoPositionWindow/M=1/R=review
DrawText 21,20,"Adjust cursor A to indicate time at which the";	DrawText 21,40," response has decayed by 90%."
Button button0,pos={80,58},size={92,20},title="Continue"; Button button0,proc=User_Continue
PauseForUser pause_for_user, review
// decay_time_0 = (xcsr(a)-response_max_time)*0.33
k2 = (xcsr(a)-response_max_time)*0.33

//duplicate /o w_in w_t
//w_t = x
//duplicate /o w_t fit_w_average_response
k4 = mean(w_in,0,uncage_time)
// k0 = v_amplitude_0
//k1 = uncage_time + v_delay_to_response_start
// k2 = decay_time_0
// k3 = rise_time_0
W_Coef = {k0, k1, k2, k3, k4}
end

function Uncaging_Analysis(data_wave_list)
	String data_wave_list
	String uncaging_response_wave_name	//name of the uncaging response wave
	Prompt uncaging_response_wave_name,"uncaging response wave name",popup,data_wave_list+"Some other wave..."
	String uncaging_power_wave_name	//name of the uncaging power wave
	Prompt uncaging_power_wave_name,"uncaging power wave name",popup,data_wave_list+"Some other wave..."
	Variable n_uncaging_pulses //variable that holds the number of uncaging pulses found in uncaging_power_wave
	Variable fit_range = 0.055	//length of the time window used for fitting double exponential
	Prompt fit_range,"Size of time window for fitting"
	Variable i = 0, j = 0	//dummy variable for iteration control
	Variable decay_time_0 = 0.008	//initial estimate for uncaging response decay time
	Prompt decay_time_0,"response decay time initial estimate"
	Variable rise_time_0 = 0.001	//initial estimate for uncaging response rise time
	Prompt rise_time_0,"response decay time initial estimate"
	Variable v_delay_to_response_start = 0
	Variable v_amplitude_0
	Variable v_amplitude_0_window = 0.001
	Prompt v_amplitude_0_window,"Time windown size for averaging to estimate peak amplitude initial estimate"
	Variable t0_box_constraint = 0.005 //size of box constraint on t0 parameter
	// Variable delay_time_to_max_response = 0.007	//not used
	Variable y0_time_window = 0.01	//time window before uncaging pulse used to estimate y0
	// Variable peak_range = 0.002	//not used
	Variable uncage_time	//time of uncaging event for current fit
	// Variable response_max_amplitude	//not used
	Variable response_max_time	//time of peak amplitude
	Variable user_response	//not used
	// Variable V_FitMaxIters = 100	//not used
	Variable threshold	//threshold value used while examining uncaging_power_wave to determine whether an algorithmically found peak is a true pulse
	// Variable V_FitError	//not used
	// Variable cursor_a, cursor_b	//not used
	Make /O/D w_coef = NaN	//not used
	Make /O/D w_sigma = NaN	//not used
	Variable V_LevelX = 0	//return value from Igor's built-in findlevel function, used for finding uncaging pulses
	// Variable y_range	//not used
	Variable last_uncage_time	//variable to temporarily hold last uncaging time
	// Variable fit_action	//not used
	// Variable y_min	//not used
	// Variable y_max	//not used
	// Variable peak_time	//not used
	Variable k0, k1, k2, k3, k4	//initial estimates for fit parameters
	Variable k0_ub // upper bounds for fit parameters
	Variable k0_lb // lower bounds for fit parameters
	Variable fit_start, fit_stop	//variables defining start and stop of the fit window
	Variable V_AbortCode = 0, V_FitError=0 //Igor environment variable indicating fitting abort condition
	Variable v_n_fit_points = Inf// number of points included in fit
	Variable n_false_replicates = 50 //number of "fake responses" to fit, per stimulus
	Variable inter_fit_time //time between end of a fit and start of the next fit for consecutive fits to "real" stimuli
	Variable inter_false_fit_time //time between fake responses
	Variable response_max_time_0, v_delay_to_response_start_0 //intial parameter estimates
	Variable v_amplitude


	String wave_data_list =  "w_amplitude;w_t0;w_decay_time;w_rise_time;w_uncage_time;w_onset_delay;w_y0;w_amplitude_se;"
	Variable n_data_waves = itemsinlist(wave_data_list)
	doPrompt "",uncaging_response_wave_name,uncaging_power_wave_name
	if(!(cmpstr(uncaging_response_wave_name, "Some other wave..." )*cmpstr(uncaging_power_wave_name, "Some other wave..." )))
	Prompt uncaging_response_wave_name,"uncaging response wave name",popup,wavelist("*",";","")
	Prompt uncaging_power_wave_name,"uncaging power wave name",popup,wavelist("*",";","")
	doPrompt "",uncaging_response_wave_name,uncaging_power_wave_name
	endif

	wave uncaging_response_wave = $uncaging_response_wave_name
	wave uncaging_power_wave = $uncaging_power_wave_name

  duplicate /o uncaging_response_wave w_uncage_response
duplicate /o uncaging_power_wave w_uncage_power

	//before we can make waves to put the fit parameters, we need to know how long to make them
	//loop through uncaging events to count the number of uncaging pulses
	//we find the first two pulses manually before entering the do loop
	//as we treat the first pulse slightly different than the rest, treating the first pulses outside of the loop allows us to avoid using an if-then construct inside the loop
	threshold = 0.8*wavemax(uncaging_power_wave)
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
	make /o/n=(n_uncaging_pulses) w_amplitude_0
	make /o/n=(n_uncaging_pulses) w_amplitude_0_alt
	make /o/n=(n_uncaging_pulses) w_amplitude_0_se
//find number of points in the uncaging response
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
//		print	x2pnt(uncaging_response_wave, fit_start ), x2pnt(uncaging_response_wave, fit_stop ), (x2pnt(uncaging_response_wave, fit_stop )-x2pnt(uncaging_response_wave, fit_start ))
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

// to get user input for initial paramters we call User_Define_Initial_Estimates
// the function has inputs (w_in,w_coef,uncage_time,y0_time_window,v_amplitude_0_window)
// //because we have aligned the resonse traces at the begining of the fit window, the uncaging pulse occurs at time = y0_time_window

User_Define_Initial_Estimates(w_avg_response,w_coef,y0_time_window,y0_time_window,v_amplitude_0_window, response_max_time, v_delay_to_response_start)
response_max_time_0 = response_max_time
v_delay_to_response_start_0 = v_delay_to_response_start



V_AbortCode = 0
V_FitError = 0
Make/O/T/N=2 T_Constraints
T_Constraints[0] = {"K0 > -30","K0 < 30","K1 > "+num2str(y0_time_window-t0_box_constraint),("K1 < "+num2str(y0_time_window+t0_box_constraint)),"K2 > 0.0001","K2 < 0.1","K3 > 0.0001","K3 < 0.1","K4 > -100","K4 < 100"}
FuncFit/N/Q/H="00000" /NTHR=0 DiffTwoExp2 W_coef  w_avg_response /D /C=T_Constraints
v_delay_to_response_start_0 = w_coef[1] - y0_time_window
v_amplitude_0 = w_coef[0]
decay_time_0 = w_coef[2]
rise_time_0 = w_coef[3]

T_Constraints[0] = {"K0 > -30","K0 < 30","K1 > "+num2str(y0_time_window),("K1 < "+num2str(y0_time_window+t0_box_constraint)),"K2 > 0.001","K2 < 0.1","K3 > 0.001","K3 < 0.1","K4 > -100","K4 < 100"}
	// do fits
	Prompt user_response, "Is this fit good?", popup, "Yes: Save fit;No: Do a Refit; No response: Save zero; Too noisy, save NaN"



	user_response = 2
	for(i=0;i < n_uncaging_pulses;i+=1)	// for1

		//uncage_time = fit_start + y0_time_window
		uncage_time = w_uncage_time[i]
		fit_start = uncage_time - y0_time_window
		fit_stop = fit_start + fit_range
		k4 = mean(uncaging_response_wave,fit_start,uncage_time)
		// k0 = v_amplitude_0
		k0 = mean(uncaging_response_wave,(fit_start + response_max_time_0-v_amplitude_0_window),(fit_start + response_max_time_0 + v_amplitude_0_window)) - k4
		k1 = uncage_time + v_delay_to_response_start_0
		k2 = decay_time_0
		k3 = rise_time_0
		W_Coef = {k0, k1, k2, k3, k4}
		V_AbortCode = 0
		V_FitError = 0
//		T_Constraints[0] = {"K0 > -30","K0 < 30","K1 > "+num2str(uncage_time-t0_box_constraint),("K1 < "+num2str(uncage_time+t0_box_constraint)),"K2 > 0.0005","K2 < 0.05","K3 > 0.0005","K3 < 0.05","K4 > -100","K4 < 100"}

//set w_coef[1] to uncaging point corrected for shifted time window
		w_coef[1] = y0_time_window + v_delay_to_response_start_0
		FuncFit/N/Q/H="01110" /NTHR=0 DiffTwoExp2 W_coef  w2d_responses[i][] /D

		// FuncFit/N/Q/H="01110" /NTHR=0 DiffTwoExp2 W_coef  uncaging_response_wave(fit_start, fit_stop) /D
		w_amplitude_0[i] = w_coef[0]
		w_amplitude_0_se[i] = w_sigma[0]
		make /o /n=(v_n_fit_points) w_temp
		w_temp = w2d_responses[i][p]
		setscale /p x,0, dimdelta(uncaging_response_wave,0), w_temp
		w_amplitude_0_alt[i] = mean(w_temp,(response_max_time_0-v_amplitude_0_window),(response_max_time_0+v_amplitude_0_window)) - mean(w_temp,0,y0_time_window)


do
dowindow /k review
display /n=review w2d_responses[i][]
SetDrawEnv xcoord= bottom;SetDrawEnv dash= 3;DelayUpdate
DrawLine uncage_time,0,uncage_time,1
// print w_coef
FuncFit/N/Q/H="00000" /NTHR=0 DiffTwoExp2 W_coef  w2d_responses[i][] /D /C=T_Constraints
// FuncFit/N/Q/H="00000" /NTHR=0 DiffTwoExp2 W_coef  w2d_responses[i][] /D
DoUpdate
// FuncFit/N/Q/H="00000" /NTHR=0 DiffTwoExp2 W_coef  uncaging_response_wave(fit_start, fit_stop) /D /C=T_Constraints
DoPrompt "Goodness of Fit", user_response
switch(user_response)	// numeric switch
case 1:		// execute if case matches expression
	break						// exit from switch
case 2:
	// print w_coef
	User_Define_Initial_Estimates(w_temp,w_coef,y0_time_window,y0_time_window,v_amplitude_0_window, response_max_time, v_delay_to_response_start)
	// print w_coef
	break
case 3:
	// print "save pars to zero"
	w_coef = 0
	w_sigma = NaN
	break
case 4:
	// print "save pars to NaN"
	w_coef = NaN
	w_sigma = NaN
endswitch

while(user_response == 2)


		w_amplitude[i] = w_coef[0]
		w_amplitude_se[i] = w_sigma[0]
		w_t0[i] = w_coef[1]
		w_decay_time[i] = w_coef[2]
		w_rise_time[i] = w_coef[3]
		w_y0[i] = w_coef[4]
		w_onset_delay[i] = w_coef[1]-y0_time_window

		duplicate /o /r=(fit_start, fit_stop) uncaging_response_wave w_t
		w_t = x
		duplicate /o w_t w_fit
		w_fit = DiffTwoExp2(w_coef, w_t)
		w2d_fits[i][] = w_fit[+q]

	endfor												// for1



make /o/n=( (n_uncaging_pulses-1)*n_false_replicates,6) w2d_fake_pars
for(i=0;i < (n_uncaging_pulses-1);i+=1)	// for1
inter_fit_time = w_fit_start_time[i+1] - w_fit_stop_time[i]
inter_false_fit_time = inter_fit_time/n_false_replicates
for(j=0;j < n_false_replicates;j+=1)	// for2
fit_start = w_fit_stop_time[i] + j*inter_false_fit_time
uncage_time = fit_start + y0_time_window
fit_stop = fit_start + fit_range
k4 = mean(uncaging_response_wave,fit_start,uncage_time)
//k0 = v_amplitude_0*sign(gnoise(1))
k0 = mean(uncaging_response_wave,(fit_start + response_max_time_0 - v_amplitude_0_window),(fit_start + response_max_time_0 - v_amplitude_0_window)) - k4
v_amplitude = k0
k1 = uncage_time + v_delay_to_response_start
k2 = decay_time_0
k3 = rise_time_0
W_Coef = {k0, k1, k2, k3, k4}
V_AbortCode = 0
V_FitError = 0
T_Constraints[0] = {"K0 > -30","K0 < 30","K1 > "+num2str(uncage_time),("K1 < "+num2str(uncage_time+t0_box_constraint)),"K2 > 0.001","K2 < 0.1","K3 > 0.001","K3 < 0.01","K4 > -100","K4 < 100"}
//FuncFit/N/Q/H="01110" /NTHR=0 DiffTwoExp2 W_coef  uncaging_response_wave(fit_start, fit_stop) /D /C=T_Constraints
//w_amplitude_0[i] = w_coef[0]
//print w_coef
// FuncFit/N/Q/H="00000" /NTHR=0 DiffTwoExp2 W_coef  uncaging_response_wave(fit_start, fit_stop) /D /C=T_Constraints


//FuncFit/N/Q/H="01110" /NTHR=0 DiffTwoExp2 W_coef  w2d_responses[i][] /D

FuncFit/N/Q/W=2/H="01110" /NTHR=0 DiffTwoExp2 W_coef  uncaging_response_wave(fit_start, fit_stop) /D

//print w_coef
// w2d_fake_pars[i*n_false_replicates+j][] = w_coef[q]
// print w_coef
w2d_fake_pars[i*n_false_replicates+j][0,4] = w_coef[q]
w2d_fake_pars[i*n_false_replicates+j][5] =  v_amplitude

// CurveFit/NTHR=0/k={fit_start + response_max_time_0} poly_XOffset 3,  w_avg_response((fit_start + response_max_time_0 - v_amplitude_0_window),(fit_start + response_max_time_0 + v_amplitude_0_window)) /D
// w2d_fake_pars[i*n_false_replicates+j][6] = poly(w_coef,fit_start + response_max_time_0) - mean(uncaging_response_wave,fit_start,(fit_start+y0_time_window))
endfor//for2
endfor//for1

end



Function DiffTwoExp2(w,t) : FitFunc
	Wave w
	Variable t

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/    IF    ((t-t0) < 0 )
	//CurveFitDialog/        f(t) = 0
	//CurveFitDialog/    ELSE
	//CurveFitDialog/        f(t) = (gsyn*(-exp(-(t-t0)/tr)+exp(-(t-t0)/td))/(-exp(-(t0+(td*tr/(td-tr))*ln(td/tr)-t0)/tr)+exp(-(t0+(td*tr/(td-tr))*ln(td/tr)-t0)/td)) -y0)
	//CurveFitDialog/    ENDIF
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ t
	//CurveFitDialog/ Coefficients 5
	//CurveFitDialog/ w[0] = gsyn
	//CurveFitDialog/ w[1] = t0
	//CurveFitDialog/ w[2] = td
	//CurveFitDialog/ w[3] = tr
	//CurveFitDialog/ w[4] = y0
	IF    ((t-w[1]) < 0 )
		return w[4]
	ELSE
		return (w[4]+w[0]*(-exp(-(t-w[1])/w[3])+exp(-(t-w[1])/w[2]))/(-exp(-(w[1]+(w[2]*w[3]/(w[2]-w[3]))*ln(w[2]/w[3])-w[1])/w[3])+exp(-(w[1]+(w[2]*w[3]/(w[2]-w[3]))*ln(w[2]/w[3])-w[1])/w[2])))
	ENDIF
End

macro Save_Results()
setscale d,0,0,"pA",w_amplitude


Display/K=0 w_amplitude
SetAxis left *,0
ModifyGraph mode=3,marker=19;DelayUpdate
ErrorBars/T=0 w_amplitude Y,wave=(w_amplitude_se,w_amplitude_se)
make /o /n=(dimsize(w2d_fake_pars,0)) w_bs_amp0
w_bs_amp0[] = w2d_fake_pars[p][0]
sort w_bs_amp0, w_bs_amp0
setscale /i x,0,1,w_bs_amp0
setscale d,0,0,"pA",w_bs_amp0
SetDrawEnv xcoord= prel;SetDrawEnv dash= 3;DelayUpdate
SetDrawEnv ycoord= left
setdrawenv linefgc=(65535,0,0)
DrawLine 0,(w_bs_amp0(0.05)),1,(w_bs_amp0(0.05))
SetDrawEnv xcoord= prel;SetDrawEnv dash= 3;DelayUpdate
SetDrawEnv ycoord= left
setdrawenv linefgc=(0,0,0)
DrawLine 0,(w_bs_amp0(0.1)),1,(w_bs_amp0(0.1))



Display/K=0 w_amplitude_0
setscale d,0,0,"pA",w_amplitude_0
ErrorBars/T=0 w_amplitude_0 Y,wave=(w_amplitude_0_se,w_amplitude_0_se)
SetAxis left *,0
ModifyGraph mode=3,marker=19;DelayUpdate
SetDrawEnv xcoord= prel;SetDrawEnv dash= 3;DelayUpdate
SetDrawEnv ycoord= left
setdrawenv linefgc=(65535,0,0)
DrawLine 0,(w_bs_amp0(0.05)),1,(w_bs_amp0(0.05))
SetDrawEnv xcoord= prel;SetDrawEnv dash= 3;DelayUpdate
SetDrawEnv ycoord= left
setdrawenv linefgc=(0,0,0)
DrawLine 0,(w_bs_amp0(0.1)),1,(w_bs_amp0(0.1))


Display/K=0 w_amplitude_0_alt
setscale d,0,0,"pA",w_amplitude_0_alt
SetAxis left *,0
ModifyGraph mode=3,marker=19;DelayUpdate
make /o /n=(dimsize(w2d_fake_pars,0)) w_bs_amp0alt
w_bs_amp0alt[] = w2d_fake_pars[p][5]
//print mean(w_temp)
//wavestats w_bs_amp0alt
//make /o/n=1 fake_amp_avg = v_avg
//make /o/n=1 fake_amp_sd = v_sdev
SetDrawEnv xcoord= prel;SetDrawEnv dash= 3;DelayUpdate
SetDrawEnv ycoord= left
sort w_bs_amp0alt, w_bs_amp0alt
setscale d,0,0,"pA",w_bs_amp0alt
setscale /i x,0,1,w_bs_amp0alt
setdrawenv linefgc=(65535,0,0)
DrawLine 0,(w_bs_amp0alt(0.05)),1,(w_bs_amp0alt(0.05))
SetDrawEnv xcoord= prel;SetDrawEnv dash= 3;DelayUpdate
SetDrawEnv ycoord= left
setdrawenv linefgc=(0,0,0)
DrawLine 0,(w_bs_amp0alt(0.1)),1,(w_bs_amp0alt(0.1))


Make/N=100/O w_bs_amp0_Hist;DelayUpdate
Histogram/B=1 w_bs_amp0,w_bs_amp0_Hist;DelayUpdate
Display w_bs_amp0_Hist
Make/N=100/O w_bs_amp0alt_Hist;DelayUpdate
Histogram/B=1 w_bs_amp0alt,w_bs_amp0alt_Hist;DelayUpdate
appendtograph /c=(0,0,0) w_bs_amp0alt_Hist


endmacro

macro Clean_Up()
	dowindow/k graph5
	dowindow/k graph4
	dowindow/k graph3
	dowindow/k graph2
	dowindow/k graph1
	dowindow/k graph0
  dowindow /k review
	//killwindow layout0
//	killwaves /a/z
kill_wave_list("w_uncage_response;w_uncage_power;")
kill_wave_list("fit_w2d_responses;w_t;w2d_fake_pars;fit_w_response_out;w_bs_amp0;w_bs_amp0alt;w_bs_amp0_Hist;w_bs_amp0alt_Hist;")
kill_wave_list("w2d_responses;w2d_fits;w_fit;w_temp;w_avg_response;T_Constraints;fit_w_avg_response;")
kill_wave_list("w_rise_time;w_y0;w_onset_delay;w_amplitude_0;w_amplitude_0_alt;w_amplitude_0_se;")
kill_wave_list("w_fit_start_time;w_fit_stop_time;w_amplitude;w_amplitude_se;w_t0;w_decay_time;")
kill_wave_list("w_response_out;w_power_out;w_coef;W_sigma;w_uncage_time;w_fit_start_pt;")
kill_wave_list("ACH_1;ACH_3;w_uncage_time;w_refs;w_stim1;w_stim2;w_stim3;w_stim4;w_stim5;w_response_out;w_power_out;w2d_responses;w2d_stim;w_temp;w_avg_response;w_avg_power;")
	killstrings /a/z
	//killdatafolder /z root:
	// setdatafolder root:
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

macro Do_Stim_Response_Average()
	String data_wave_list = s_wavenames
	stim_response_average(data_wave_list)
	//Kill_Wave_List(data_wave_list)
endmacro




Function stim_response_average(data_wave_list)
String data_wave_list
String uncaging_response_wave_name	//name of the uncaging response wave
Prompt uncaging_response_wave_name,"uncaging response wave name",popup,data_wave_list+"Some other wave..."
String uncaging_power_wave_name	//name of the uncaging power wave
Prompt uncaging_power_wave_name,"uncaging power wave name",popup,data_wave_list+"Some other wave..."
Variable threshold, n_uncaging_pulses
	Variable V_LevelX = 0	//return value from Igor's built-in findlevel function, used for finding uncaging pulses
  Variable v_delta_t = 1E-4
Variable v_fit_range = 0.055	//length of the time window used for fitting double exponential
Variable v_n_fit_points
Variable y0_time_window = 0.01
Variable v_fit_start, v_fit_stop, v_fit_start_point, v_fit_stop_point, n_reps

doPrompt "",uncaging_response_wave_name,uncaging_power_wave_name

if(!(cmpstr(uncaging_response_wave_name, "Some other wave..." )*cmpstr(uncaging_power_wave_name, "Some other wave..." )))
Prompt uncaging_response_wave_name,"uncaging response wave name",popup,wavelist("*",";","")
Prompt uncaging_power_wave_name,"uncaging power wave name",popup,wavelist("*",";","")
doPrompt "",uncaging_response_wave_name,uncaging_power_wave_name
endif

wave uncaging_response_wave = $uncaging_response_wave_name
wave uncaging_power_wave = $uncaging_power_wave_name
setscale /p x,0,1E-4, "s", uncaging_response_wave
setscale /p x,0,1E-4, "s", uncaging_power_wave

//before we can make waves to put the fit parameters, we need to know how long to make them
//loop through uncaging events to count the number of uncaging pulses
//we find the first two pulses manually before entering the do loop
//as we treat the first pulse slightly different than the rest, treating the first pulses outside of the loop allows us to avoid using an if-then construct inside the loop
threshold = 0.8*wavemax(uncaging_power_wave)
// find rising edge in power trace greater than defined threshold
findlevel /Q/EDGE=1 /R= (V_LevelX,) uncaging_power_wave, threshold
//we find the first pulse manually before entering the do loop, so we start the counter at 1
Make/O/N=1 w_uncage_time
w_uncage_time = V_LevelX
variable i,j
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


make /wave /o /n=5 w_refs
make /o/n=5 w_stim1 = {2,3,5,8,12}
w_refs[0] = w_stim1
make /o/n=4 w_stim2 = {4,6,9,13}
w_refs[1] = w_stim2
make /o/n=3 w_stim3 = {7,10,14}
w_refs[2] = w_stim3
make /o/n=2 w_stim4 = {11,15}
w_refs[3] = w_stim4
make /o/n=2 w_stim5 = {1,16}
w_refs[4] = w_stim5


v_n_fit_points = x2pnt(uncaging_response_wave, v_fit_range)

duplicate /o /r=(0,w_uncage_time[numpnts(w_refs)]) uncaging_response_wave w_response_out
duplicate /o w_response_out w_power_out
w_power_out = 0
j = 0

for(j = 0;j<numpnts(w_refs);j+=1)//for2
n_reps = numpnts(w_refs[j])
wave w = w_refs[j]

for(i=0;i<n_reps;i+=1)//for1
make /o/n=(n_reps,v_n_fit_points) w2d_responses
make /o/n=(n_reps,v_n_fit_points) w2d_stim

v_fit_start = w_uncage_time[w[i]-1]-y0_time_window
v_fit_start_point = x2pnt(uncaging_response_wave,v_fit_start)
w2d_responses[i][] =  uncaging_response_wave[v_fit_start_point + q]
w2d_stim[i][] =  uncaging_power_wave[v_fit_start_point + q]
endfor//for1

make /o /n=(n_reps) w_temp=1
MatrixOp/O w_avg_response=w_temp^t x w2d_responses
MatrixOp/O w_avg_power=w_temp^t x w2d_stim
w_avg_response = w_avg_response/n_reps
w_avg_power = w_avg_power/n_reps
redimension /n=(v_n_fit_points) w_avg_response
redimension /n=(v_n_fit_points) w_avg_power
v_fit_start = w_uncage_time[j]-y0_time_window
v_fit_start_point = x2pnt(uncaging_response_wave,v_fit_start)
w_response_out[v_fit_start_point,v_fit_start_point+v_n_fit_points-1] = w_avg_response[p - v_fit_start_point]
w_power_out[v_fit_start_point,v_fit_start_point+v_n_fit_points-1] = w_avg_power[p - v_fit_start_point]
endfor//for2





save /T w_response_out,w_power_out as "stim_response_average.itx"

end


menu "macros"
	"Do_Uncaging_Analysis/1"
	"Save_Results/2"
	"Clean_Up/3"
	"Kill_Input_Waves"
	"Do_Stim_Response_Average"
		"Do_Save_Results2/4"
end
