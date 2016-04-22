#pragma rtGlobals=3		// Use modern global access method and strict wave access.

Function UserCursorAdjust_ContButtonProc(ctrlName) : ButtonControl
	String ctrlName
	DoWindow/K tmp_PauseforCursor
End
macro Do_Uncaging_Analysis()
String data_wave_list = s_wavenames
Uncaging_Analysis(data_wave_list)
//Kill_Wave_List(data_wave_list)
endmacro

function Uncaging_Analysis(data_wave_list)
	String data_wave_list
	 String uncaging_response_wave_name	//name of the uncaging response wave
	 Prompt uncaging_response_wave_name,"uncaging response wave name",popup,data_wave_list
	 String uncaging_power_wave_name	//name of the uncaging power wave
	 Prompt uncaging_power_wave_name,"uncaging power wave name",popup,data_wave_list
	 Variable n_uncaging_pulses //variable that holds the number of uncaging pulses found in uncaging_power_wave
	Variable fit_range = 0.055	//length of the time window used for fitting double exponential
	Prompt fit_range,"Size of time window for fitting"
	Variable i = 0	//dummy variable for iteration control
	Variable decay_time_0 = 0.008	//initial estimate for uncaging response decay time
	Prompt decay_time_0,"response decay time initial estimate"
	Variable rise_time_0 = 0.004	//initial estimate for uncaging response rise time
	Prompt rise_time_0,"response decay time initial estimate"
	Variable t0_box_constraint = 0.005 //size of box constraint on t0 parameter
	// Variable delay_time_to_max_response = 0.007	//not used
	// Variable y0_time_window = 0.01	//not used
	// Variable peak_range = 0.002	//not used
	Variable uncage_time	//time of uncaging event for current fit
	// Variable response_max_amplitude	//not used
	// Variable response_max_time	//not used
	// Variable user_response	//not used
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
	// Variable k0_0, k1_0, k2_0	//not used
	Variable fit_start, fit_stop	//variables defining start and stop of the fit window
	Variable V_AbortCode = 0 //Igor environment variable indicating fitting abort condition

	String wave_data_list =  "w_amplitude;w_t0;w_decay_time;w_rise_time;w_uncage_time;w_onset_delay;w_y0;w_amplitude_se"
	Variable n_data_waves = itemsinlist(wave_data_list)
doPrompt "",uncaging_response_wave_name,uncaging_power_wave_name
wave uncaging_response_wave = $uncaging_response_wave_name
wave uncaging_power_wave = $uncaging_power_wave_name



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

//create waves to store fitted parameters
	i=0
	do // do1
		//if (!waveexists($stringfromlist(i,wave_data_list)))//if1
			Make/O/N=(n_uncaging_pulses) $stringfromlist(i,wave_data_list)
		//endif//if1
		i = i + 1
	while(i < n_data_waves) //do1


fit_start = 0.015
uncage_time = fit_start + 0.01
fit_stop = fit_start + fit_range
duplicate/o /r=(fit_start, fit_stop) uncaging_response_wave w_response_temp
W_coef = {1, uncage_time, decay_time_0, rise_time_0, mean(w_response_temp)}
V_AbortCode = 0
Make/O/T/N=2 T_Constraints
//T_Constraints[0] = {("K1 > "+num2str(uncage_time-t0_box_constraint)),("K1 < "+num2str(uncage_time+t0_box_constraint))}
T_Constraints[0] = {"K0 > -5","K0 < 5","K1 > 0","K1 < 1","K2 > 0.0001","K2 < 0.1","K3 > 0.0001","K3 < 0.1","K4 > -100","K4 < 100"}
FuncFit/N/Q/H="00000" /NTHR=0 DiffTwoExp2 W_coef  w_response_temp /D /C=T_Constraints
//FuncFit/N/Q/H="00000" /NTHR=0 DiffTwoExp2 W_coef  w_response_temp /D

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
endmacro

macro Clean_Up()
dowindow/k graph5
dowindow/k graph4
dowindow/k graph3
dowindow/k graph2
dowindow/k graph1
dowindow/k graph0
//killwindow layout0
killwaves /a/z
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

menu "macros"
	"Do_Uncaging_Analysis/1"
	"Save_Results/2"
	"Clean_Up/3"
	"Kill_Input_Waves"
end
