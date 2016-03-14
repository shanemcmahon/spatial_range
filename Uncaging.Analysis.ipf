#pragma rtGlobals=3		// Use modern global access method and strict wave access.
save /b/t/o/p=output_path  wavelist("*_u*",";","")+wavelist("*_p*",";","")+wave_data_list as replacestring(".itx", S_fileName, "_out.itx")

//*******************************************************************************************************************************

Function PopupChoice ()

	String tracename
	Variable Choose = 4
	Prompt traceName,"Trace",popup,TraceNameList("",";",1)
	Prompt Choose, "Is this fit good?", popup, "No: Do a Refit; No response: Save zero; Too noisy, save NaN; Yes: Save fit"
	DoPrompt "Goodness of Fit", Choose
	return choose
End

function insert_fit_coefs(w_coef, amplitudeData, T0Data, TauDecayData, TauRiseData, OffsetData, TimeDiffUncgpnt_Response, Amplitude_SD, uncgpnt, w_sigma)
	wave w_coef, amplitudeData, t0data, TauDecayData, TauRiseData, OffsetData, TimeDiffUncgpnt_Response, Amplitude_SD, w_sigma
	variable uncgpnt

	// InsertPoints numpnts($stringfromlist(i,wave_data_list)), 1, $stringfromlist(i,wave_data_list)
	// $stringfromlist(i,wave_data_list)[numpnts($stringfromlist(i,wave_data_list))-1] = w_coef[i]
	// string /g wave_data_list =  "AmplitudeData;T0Data;TauDecayData;TauRiseData;OffsetData;uncgpntData;TimeDiffUncgpnt_Response;Amplitude_SD"
	InsertPoints numpnts(amplitudeData), 1, amplitudeData
	amplitudeData[numpnts(amplitudeData)-1] = W_coef[0]

	InsertPoints numpnts(T0Data), 1, T0Data
	T0Data[numpnts(T0Data)-1] = W_coef[1]

	InsertPoints numpnts(TauDecayData), 1, TauDecayData
	TauDecayData[numpnts(TauDecayData)-1] = W_coef[2]

	InsertPoints numpnts(TauRiseData), 1,TauRiseData
	TauRiseData[numpnts(TauRiseData)-1] = W_coef[3]

	InsertPoints numpnts(OffsetData), 1, OffsetData
	OffsetData[numpnts(OffsetData)-1] = W_coef[4]

	InsertPoints numpnts(TimeDiffUncgpnt_Response), 1, TimeDiffUncgpnt_Response
	TimeDiffUncgpnt_Response[numpnts(TimeDiffUncgpnt_Response)-1] = ((W_coef[1]) - uncgpnt)

	InsertPoints numpnts(Amplitude_SD), 1, Amplitude_SD
	Amplitude_SD[numpnts(Amplitude_SD)-1] = W_Sigma[0]
end

Function Refit()
	Variable better = 2
	Prompt better, "Better?", popup, "yes; no"
	String inputstring = "Improved Fit"
	DoPrompt/HELP="if you click cancel it will return yes and save" inputstring, better
	return better
End

Function UserCursorAdjust_ContButtonProc(ctrlName) : ButtonControl
	String ctrlName
	DoWindow/K tmp_PauseforCursor
End


function do_fit(traceunc, xmin, xmax, w_coef)
	wave w_coef
	string traceunc
	variable xmin, xmax
	variable v_abortcode = 0
	try
		//		Display/N=Checking $traceunc;
		//SetAxis bottom (uncgpnt-0.01),(uncgpnt+(2*fitrange))
		FuncFit/N/Q/H="00000" /NTHR=0 DiffTwoExp2 W_coef  $traceunc(xmin,xmax) /D
		//		AppendToGraph /c=(0,0,0) $("fit_"+traceunc)
		//		ModifyGraph rgb($("fit_"+traceunc))=(0,0,0)
		//		abort("1")
	catch
		W_coef = NaN
		Duplicate/O/R=(xmin,xmax) $traceunc, $("fit_"+traceunc)
	endtry

end

Macro UncagingAnalysis (traceunc, uncageinterval, fitrange)
	String traceunc = "ach_1"
	Variable uncageinterval=1, fitrange=0.04
	// Prompt traceunc, "Response Trace", popup, WaveList("*", ";", "")
	// Prompt uncageinterval, "Interval between uncaging points (s)"
	// Prompt fitrange, "Fitting range (s)"
	String/G tracepower = "ach_3"
	Variable i=0
	Variable Td = 0.004							//Tau decay estimate
	Variable Tr = 0.002						//Tau rise estimate
	Variable peak_loc = .007
	Variable offsetrange = 0.01					//10ms averaged for offset of trace
	Variable peakrange = 0.002				//1ms averaged for offset of peak
	Variable/G uncgpnt						//uncaging point and it's end
	Variable maximum							//maximum of each uncaging power pulse from trace
	Variable/G betterreturn
	Variable V_FitMaxIters = 100			//Number of iterations FitFunc can go through
	Variable/G Threshold
	Variable/G V_FitError					//If there's a problem e.g. Singular Matrix Error, then saves specifies action (save NaN)
	Variable/G cursorA, cursorB
	Make/O/D W_coef = NaN					//For holding DiffTwoExp coefficients
	Make/O/D W_sigma = NaN					//For holding DiffTwoExp coefficients
	variable v_levelx = 0
	variable y_range = 0
	variable last_uncage_time = 0
	variable t_max_response
	variable fit_action
	variable y_min
	variable y_max
  variable peak_time, k0_0, k1_0, k2_0,fit_start, fit_stop
	string /g wave_data_list =  "AmplitudeData;T0Data;TauDecayData;TauRiseData;uncgpntData;TimeDiffUncgpnt_Response;OffsetData;Amplitude_SD"
	variable n_data_waves = itemsinlist(wave_data_list)

	i=0
	do // do1
		if (!waveexists($stringfromlist(i,wave_data_list)))
			Make/O/N=0 $stringfromlist(i,wave_data_list)
		endif
		i = i + 1
	while(i < n_data_waves) //do1



	//--------------------------------Start loop
	// Initial loop through the data to find an appropriate y-range for displaying plots
	i = 0
	Duplicate/O $traceunc, smoothed
	Smooth/B=7 1, smoothed
	//use find level to find uncaging pulses in tracepower
	// set threshold value to look for peaks so that findlevel does not find false peaks in the noise
	threshold = 0.8*wavemax($tracepower)
	findlevel /Q /EDGE=1 /R= (v_levelx,) $tracepower, threshold
  uncgpnt = v_levelx
  String/G fitwave = nameofwave($traceunc) + "_u" + num2str(i+1)
  String/G tracecopy = nameofwave($traceunc) + "_P" + num2str(i+1)
  string/g uncagecopy = nameofwave($traceunc) + "_pockel" + num2str(i+1)
  duplicate/o /r=((uncgpnt-fitrange/3), (uncgpnt+fitrange)) $tracepower $uncagecopy
  duplicate/o /r=((uncgpnt-fitrange/3), (uncgpnt+fitrange)) $traceunc $tracecopy
  duplicate /o $tracecopy average_trace
	findlevel /Q/EDGE=2 /R= (v_levelx,) $tracepower, threshold
	findlevel /Q/EDGE=1 /R= (v_levelx,) $tracepower, threshold
i=1
  uncgpnt = v_levelx
  String/G fitwave = nameofwave($traceunc) + "_u" + num2str(i+1)
  String/G tracecopy = nameofwave($traceunc) + "_P" + num2str(i+1)
  string/g uncagecopy = nameofwave($traceunc) + "_pockel" + num2str(i+1)
  duplicate/o /r=((uncgpnt-fitrange/3), (uncgpnt+fitrange)) $tracepower $uncagecopy
  duplicate/o /r=((uncgpnt-fitrange/3), (uncgpnt+fitrange)) $traceunc $tracecopy
  average_trace = average_trace + $tracecopy
	last_uncage_time = v_levelx
	// v_levelx is now set to the rising phase of the second uncaging event
	// the interval (0,v_levelx) contains the first response
	// y_range = wavemax(smoothed,0,v_levelx)-wavemin(smoothed,0,v_levelx)
	y_range = wavemax(smoothed,(last_uncage_time-0.5*fitrange),(last_uncage_time+fitrange))-wavemin(smoothed,(last_uncage_time-0.5*fitrange),(last_uncage_time+fitrange))
	do //do2
		// find the next falling phase followed bz the next rising phase
		findlevel /Q/EDGE=2 /R= (v_levelx,) $tracepower, threshold
		findlevel /Q/EDGE=1 /R= (v_levelx,) $tracepower, (threshold)
		// the interval (last_uncage_time,v_levelx) now contains the i+1th uncaging response
		// set y_range to the larger of the previous y_range or the difference between the max/min during the current uncaging event time window
		// y_range = max(y_range,wavemax(smoothed,last_uncage_time,v_levelx)-wavemin(smoothed,last_uncage_time,v_levelx))
		y_range = max(y_range,(wavemax(smoothed,(last_uncage_time-0.5*fitrange),(last_uncage_time+fitrange))-wavemin(smoothed,(last_uncage_time-0.5*fitrange),(last_uncage_time+fitrange))))
		// when the last uncaging pulse has passed, findlevel will fail to find a rising edge on the inverval (v_levelx,inf)
		// at this point, findlevel returns NaN which signals that the last uncaging event has been processed
		// NaN==NaN evaluates to logically false, therefore, NaN is detected by Igor's built in numtype function, which returns 0 for a numeric value, 1 for inf or 2 for NaN, therefore numtype(v_levelx) evaluates to logically true when v_levelx == NaN
		if(numtype(v_levelx))
			break
		endif
		last_uncage_time = v_levelx
		i = i + 1
    uncgpnt = v_levelx
    String/G fitwave = nameofwave($traceunc) + "_u" + num2str(i+1)
    String/G tracecopy = nameofwave($traceunc) + "_P" + num2str(i+1)
    string/g uncagecopy = nameofwave($traceunc) + "_pockel" + num2str(i+1)
    duplicate/o /r=((uncgpnt-fitrange/3), (uncgpnt+fitrange)) $tracepower $uncagecopy
    duplicate/o /r=((uncgpnt-fitrange/3), (uncgpnt+fitrange)) $traceunc $tracecopy
    average_trace = average_trace + $tracecopy
	while(1)//do2
	i = 0; v_levelx = 0
	do //do3
		V_FitError = 0;

		// find uncaging event from pockels cell signal
		findlevel /Q/EDGE=1 /R= (v_levelx,) $tracepower, threshold
		// when the last uncaging pulse has passed, findlevel will fail to find a rising edge on the inverval (v_levelx,inf)
		// at this point, findlevel returns NaN which signals that the last uncaging event has been processed
		// NaN==NaN evaluates to logically false, therefore, NaN is detected by Igor's built in numtype function, which returns 0 for a numeric value, 1 for inf or 2 for NaN, therefore numtype(v_levelx) evaluates to logically true when v_levelx == NaN

		if(numtype(v_levelx))
			break
		endif
		uncgpnt = v_levelx
		findlevel /Q/EDGE=2 /R= (v_levelx,) $tracepower, threshold
		InsertPoints numpnts (uncgpntData), 1, uncgpntData
		uncgpntData[numpnts(uncgpntData)] = uncgpnt

		// uncgpnt = v_levelx
    String/G fitwave = nameofwave($traceunc) + "_u" + num2str(i+1)
    String/G tracecopy = nameofwave($traceunc) + "_P" + num2str(i+1)
    string/g uncagecopy = nameofwave($traceunc) + "_pockel" + num2str(i+1)
    // duplicate/o /r=((uncgpnt-fitrange/3), (uncgpnt+fitrange)) $tracepower $uncagecopy
    // duplicate/o /r=((uncgpnt-fitrange/3), (uncgpnt+fitrange)) $traceunc $tracecopy

if(i==0)
// display average_trace
K0 =  mean(average_trace,(uncgpnt - offsetrange),uncgpnt)
Loess/pass=1 /DEST=temp srcWave=average_trace
Display/N=checkavg average_trace;
SetDrawEnv xcoord= bottom;SetDrawEnv dash= 3;DelayUpdate
DrawLine uncgpnt,0,uncgpnt,1
// AppendToGraph /c=(0,0,0) $("fit_"+traceunc)
// AppendToGraph /C = (0,0,0)$fitwave;
//SetAxis/W=checkavg left (wavemax(average_trace)-y_range), (wavemax(average_trace))
NewPanel/K=2 as "Pause for cursor";	DoWindow/C tmp_PauseforCursor; AutoPositionWindow/M=1/R=checkavg
DrawText 21,20,"Adjust cursor A to estimate y0 and then";	DrawText 21,40,"Click Continue."
Button button0,pos={80,58},size={92,20},title="Continue"; Button button0,proc=UserCursorAdjust_ContButtonProc
ShowInfo/CP=0/W=checkavg; PauseForUser tmp_PauseforCursor, checkavg
Wavestats/Q/R=(xcsr(a), xcsr(b)) average_trace;	w_coef[4]=V_avg
NewPanel/K=2 as "Pause for cursor"; DoWindow/C tmp_PauseforCursor; AutoPositionWindow/M=1/R=checkavg

DrawText 21,20,"Adjust cursor A to estimate response"; DrawText 21,40,"Click Continue."; Button button0,pos={80,58},size={92,20},title="Continue"
Button button0,proc=UserCursorAdjust_ContButtonProc; PauseForUser tmp_PauseforCursor, checkavg

Wavestats/Q/R=(xcsr(A), xcsr(B)) average_trace;	w_coef[0] = v_avg - w_coef[4];	w_coef[1] = uncgpnt; t_max_response = (xcsr(a)+xcsr(b))/2; w_coef[3] = (t_max_response - uncgpnt)/3

NewPanel/K=2 as "Pause for cursor";	DoWindow/C tmp_PauseforCursor;	AutoPositionWindow/M=1/R=checkavg
DrawText 21,20,"Adjust cursor A to estimate decay time";	DrawText 21,40,"Click Continue.";	Button button0,pos={80,58},size={92,20},title="Continue"
Button button0,proc=UserCursorAdjust_ContButtonProc;	PauseForUser tmp_PauseforCursor, checkavg
w_coef[2] = (xcsr(a) - t_max_response)/2

NewPanel/K=2 as "Pause for cursor";	DoWindow/C tmp_PauseforCursor;	AutoPositionWindow/M=1/R=checkavg
DrawText 21,20,"Adjust cursors to define fit window";	DrawText 21,40,"Click Continue."; Button button0,pos={80,58},size={92,20},title="Continue"
Button button0,proc=UserCursorAdjust_ContButtonProc;	PauseForUser tmp_PauseforCursor, checkavg
cursorA = xcsr(A);	cursorB = xcsr(B)
fit_start = xcsr(A)-uncgpnt;	fit_stop = xcsr(B)-uncgpnt

DoWindow/K checkavg
temp = (abs(temp - K0))
wavestats /q temp
peak_time = v_maxloc
K1 = temp(v_maxloc)-k0
K2 = td
k0_0 = k0
k1_0 = k1
k2_0 = k2
// obtain initial estimates for amplitude, y0, and decay period from a fit to a single exponential
//CurveFit/N/Q/NTHR=0 exp_XOffset average_trace(uncgpnt + peak_loc,(uncgpnt + peak_loc + fitrange)) /D
CurveFit/N/Q/NTHR=0 exp_XOffset average_trace(t_max_response,fit_stop+uncgpnt) /D

// if estimated decay from single exponential is very large, some numerical instability may result
// typical decay times are on the order of ms, if the estimated decay is >0.5, then replace by the arbitrary cutoff value .02 to improve stability
if(W_coef[2]>.5)
  W_coef[2]=.02
endif
// The decay time should be positive, the rough estimation procedure may occasionaly produce a negative estimate, particularly when there is no response
// A negative value for the rise time causes numerical errors in DiffTwoExp2, so we will replace it with an arbitrary value when it occurs
if(W_coef[2]<0)
  W_coef[2] = 0.02
endif
if(v_fiterror)
  // if the exponential fit causes an error, then use a generic estimate for initial parameters
  v_fiterror = 0
  // W_coef = {mean($traceunc,(uncgpnt - offsetrange),uncgpnt),-1E-12,.04}
      W_coef = {k1_0, uncgpnt, k2_0, ((peak_time-uncgpnt)*0.33), k0_0}
endif
W_coef = {(W_coef[1]), uncgpnt, W_coef[2], ((peak_time-uncgpnt)*0.33), W_coef[0]}
	// do_fit("average_trace", (uncgpnt-fitrange/3), (uncgpnt+fitrange), w_coef)
  	do_fit("average_trace", fit_start+uncgpnt, fit_stop+uncgpnt, w_coef)
  duplicate /o w_coef w_coef_0
endif

    K0 =  mean($traceunc,(uncgpnt - offsetrange),uncgpnt)
		Loess/pass=1 /DEST=temp srcWave=$tracecopy
		temp = (abs(temp - K0))
    wavestats /q temp
    peak_time = v_maxloc
    K1 = temp(v_maxloc)-k0
    K2 = k2_0
    W_coef = {k1, uncgpnt, W_coef_0[2], ((peak_time-uncgpnt)*0.33), k0}
      // W_coef = {k1_0, uncgpnt, k2_0, ((peak_time-uncgpnt)*0.33), k0_0}
    // k0_0 = k0
    // k1_0 = k1
    // k2_0 = k2
    // obtain initial estimates for amplitude, y0, and decay period from a fit to a single exponential
    //CurveFit/N/Q/NTHR=0 exp_XOffset $traceunc(uncgpnt + peak_loc,(uncgpnt + peak_loc + fitrange)) /D
        CurveFit/N/Q/NTHR=0 exp_XOffset $traceunc(peak_time, fit_stop+uncgpnt) /D
    // if estimated decay from single exponential is very large, some numerical instability may result
    // typical decay times are on the order of ms, if the estimated decay is >0.5, then replace by the arbitrary cutoff value .02 to improve stability
    if(W_coef[2]>.5)
      W_coef[2]=.02
    endif
    // The decay time should be positive, the rough estimation procedure may occasionaly produce a negative estimate, particularly when there is no response
    // A negative value for the rise time causes numerical errors in DiffTwoExp2, so we will replace it with an arbitrary value when it occurs
    if(W_coef[2]<0)
      W_coef[2] = 0.02
    endif
    if(v_fiterror)
      // if the exponential fit causes an error, then use a generic estimate for initial parameters
      v_fiterror = 0
      // W_coef = {mean($traceunc,(uncgpnt - offsetrange),uncgpnt),-1E-12,.04}
          //W_coef = {k1_0, uncgpnt, k2_0, ((peak_time-uncgpnt)*0.33), k0_0}
          w_coef = w_coef_0
    endif
		W_coef = {(W_coef[1]), uncgpnt, W_coef[2], ((peak_time-uncgpnt)*0.33), W_coef[0]}
    		// W_coef = {(W_coef[1]), uncgpnt, W_coef[2], w_coef_0[3], W_coef[0]}
    // W_coef = {k1, uncgpnt, k2, ((peak_time-uncgpnt)*0.33), k0}



		//-------------------------------------Name graphs and Display



//    Display/N=Checking $traceunc;
    Display/N=Checking $tracecopy;
		y_min = (wavemax($traceunc,(uncgpnt-0.01),(uncgpnt+(2*fitrange)))-y_range)
		y_max = wavemax($traceunc,(uncgpnt-0.01),(uncgpnt+(2*fitrange)))
		SetAxis/W=Checking left y_min, y_max;
		SetAxis bottom fit_start+uncgpnt,fit_stop+uncgpnt
		SetDrawEnv xcoord= bottom;SetDrawEnv dash= 3;DelayUpdate
		DrawLine uncgpnt,0,uncgpnt,1
    // fit difference in exponential function
		// do_fit(traceunc, (uncgpnt-fitrange/3), (uncgpnt+fitrange), w_coef)
    do_fit(traceunc, fit_start+uncgpnt, fit_stop+uncgpnt, w_coef)
		print w_coef[0], w_sigma[0]

		duplicate /o fit_ach_1 $fitwave
		AppendToGraph /c=(0,0,0) $("fit_"+traceunc)
   appendtograph /r /t /c=(0,0,65535) fit_average_trace
		//---------------------------------------------Start of if-loop popup menu: whether to save data=good fits or NaN=bad fits
		//---------------------------------------------
		fit_action = PopupChoice ()
		// if fit_action == 4, then the fit is good and coefficients are saved, otherwise fit is bad, we save, 0, NaN or do a refit
		if (fit_action == 2)		//Nothing there: save as zero for no response
			w_coef = 0
			w_sigma = NaN
		endif
		if (fit_action == 3)		//Nothing there: noise! Save NaN
			w_coef = NaN
			w_sigma = NaN
		endif
		if(fit_action != 1)
			insert_fit_coefs(w_coef, amplitudeData, T0Data, TauDecayData, TauRiseData, OffsetData, TimeDiffUncgpnt_Response, Amplitude_SD, uncgpnt, w_sigma)
		endif
		if(fit_action == 1)
			//No, not good fit. Get fit points from cursors
			DoWindow/K Checking
			Duplicate/O/R=((uncgpnt-(0.5*fitrange)),(uncgpnt+(2*fitrange))) $traceunc, $tracecopy
			// Duplicate/O/R=(uncgpnt,(uncgpnt+fitrange)) fit_wave1, $fitwave
			loess /N=13 /ord=2 /dest=$tracecopy srcwave=$tracecopy
			Display/N=Checking $tracecopy;
			SetDrawEnv xcoord= bottom;SetDrawEnv dash= 3;DelayUpdate
			DrawLine uncgpnt,0,uncgpnt,1
			// AppendToGraph /c=(0,0,0) $("fit_"+traceunc)
			// AppendToGraph /C = (0,0,0)$fitwave;
			SetAxis/W=Checking left (wavemax($tracecopy)-y_range), (wavemax($tracecopy))
			NewPanel/K=2 as "Pause for cursor";	DoWindow/C tmp_PauseforCursor; AutoPositionWindow/M=1/R=Checking
			DrawText 21,20,"Adjust cursor A to estimate y0 and then";	DrawText 21,40,"Click Continue."
			Button button0,pos={80,58},size={92,20},title="Continue"; Button button0,proc=UserCursorAdjust_ContButtonProc
			ShowInfo/CP=0/W=Checking; PauseForUser tmp_PauseforCursor, Checking
			Wavestats/Q/R=(xcsr(a), xcsr(b)) $traceunc;	w_coef[4]=V_avg
			NewPanel/K=2 as "Pause for cursor"; DoWindow/C tmp_PauseforCursor; AutoPositionWindow/M=1/R=Checking

			DrawText 21,20,"Adjust cursor A to estimate response"; DrawText 21,40,"Click Continue."; Button button0,pos={80,58},size={92,20},title="Continue"
			Button button0,proc=UserCursorAdjust_ContButtonProc; PauseForUser tmp_PauseforCursor, Checking

			Wavestats/Q/R=(xcsr(A), xcsr(B)) $traceunc;	w_coef[0] = v_avg - w_coef[4];	w_coef[1] = uncgpnt; t_max_response = (xcsr(a)+xcsr(b))/2; w_coef[3] = (t_max_response - uncgpnt)/3

			NewPanel/K=2 as "Pause for cursor";	DoWindow/C tmp_PauseforCursor;	AutoPositionWindow/M=1/R=Checking
			DrawText 21,20,"Adjust cursor A to estimate decay time";	DrawText 21,40,"Click Continue.";	Button button0,pos={80,58},size={92,20},title="Continue"
			Button button0,proc=UserCursorAdjust_ContButtonProc;	PauseForUser tmp_PauseforCursor, Checking
			w_coef[2] = (xcsr(a) - t_max_response)/2

			NewPanel/K=2 as "Pause for cursor";	DoWindow/C tmp_PauseforCursor;	AutoPositionWindow/M=1/R=Checking
			DrawText 21,20,"Adjust cursors to define fit window";	DrawText 21,40,"Click Continue."; Button button0,pos={80,58},size={92,20},title="Continue"
			Button button0,proc=UserCursorAdjust_ContButtonProc;	PauseForUser tmp_PauseforCursor, Checking
			cursorA = xcsr(A);	cursorB = xcsr(B)

			DoWindow/K Checking

			//---------------------------------------------Display and check if good fit
			String/G fitwave = nameofwave($traceunc) + "_u" + num2str(i+1)
			String/G tracecopy = nameofwave($traceunc) + "_P" + num2str(i+1)

			Display/N=Checking $tracecopy;

			SetAxis/W=Checking left (wavemax($traceunc,cursorA,cursorB)-y_range), wavemax($traceunc,cursorA,CursorB);
			SetAxis bottom CursorA,CursorB
			do_fit(traceunc, CursorA, CursorB, w_coef)
			duplicate/o /r=((uncgpnt-fitrange/3), (uncgpnt+fitrange)) $traceunc $tracecopy
			duplicate/o fit_ach_1 $fitwave
			string/g uncagecopy = nameofwave($traceunc) + "_pockel" + num2str(i+1)
			duplicate/o /r=((uncgpnt-fitrange/3), (uncgpnt+fitrange)) $tracepower $uncagecopy
			AppendToGraph /c=(0,0,0) $("fit_"+traceunc)
			print w_coef[0], w_sigma[0]
			betterreturn = Refit()



			if (betterreturn == 2)//-------------------------------------If fit not good, saves in waves as NaN
				w_coef = NaN
				w_sigma = NaN
			endif
			insert_fit_coefs(w_coef, amplitudeData, T0Data, TauDecayData, TauRiseData, OffsetData, TimeDiffUncgpnt_Response, Amplitude_SD, uncgpnt, w_sigma)
		endif //end if fit_action == 1

		DoWindow/K Checking
		i += 1
	while(1) //do3
	// Save/T AmplitudeData, T0Data, TauDecayData, TauRiseData, TimeDiffUncgpnt_Response, OffsetData, Amplitude_SD as replacestring(".itx", S_fileName, "_out.itx")
	newpath /o output_path s_path
//saveexperiment /p = output_path as  s_filename[0,(strsearch(s_filename,".",0)-1)]+".pxp"
//save /b/t/o/p=output_path  wavelist("*_u*",";","")+wavelist("*_p*",";","")+wave_data_list+";average_trace;fit_average_trace" as replacestring(".itx", S_fileName, "_out.itx")
	NewPath/o referencepath  s_path+"References"
plot_waves()

	// KillWaves W_coef, smoothed, fit
	// KillVariables/A;	KillStrings/A;	KillWaves /A
EndMacro

macro plot_waves()
//	newpath /o output_path s_path
//	NewPath/o referencepath  s_path+"References"
ImageLoad/P=referencepath/T=tiff Indexedfile(referencepath, 0,".tif")
//	duplicate /o AmplitudeData tmp
//	duplicate /o Amplitude_SD tmp2
//	tmp = AmplitudeData[numpnts(amplitudedata)-1-p]
//	tmp2 = Amplitude_sd[numpnts(amplitudedata)-1-p]
//	amplitudedata = -tmp
//	Amplitude_sd = tmp2
if(!strlen(winlist("layout0",";","")))
	newLayout /P=Landscape/w=(0,0,1000,650)
	endif
	modifylayout mag=1
	NewImage/K=0  $StringFromList(0,WaveList("Trigger*",";",""))
	ModifyGraph width=210,height=210
	ModifyGraph width=210,height={Aspect,1}
	AppendToLayout/T Graph0
	ModifyLayout left(Graph0)=0,top(Graph0)=0
	display ach_1
	Label left "I (pA)"
	ModifyGraph lsize=0.1
	ModifyGraph width=500,height=167.5
	AppendToLayout/T Graph1
	ModifyLayout left(Graph1)=230,top(Graph1)=20
	display ach_1_p10
	ModifyGraph lsize=0.1
	AppendToGraph /C = (0,0,0) ach_1_u10
	appendtograph /r /c=(0,0,65535) ACH_1_pockel10
	Label left "I (pA)"
	Label right "V\\BPockels \\M(V)"
	ModifyGraph width=100,height=125
	AppendToLayout/T Graph2
	ModifyLayout left(Graph2)=35,top(Graph2)=240

	display ach_1_p5
	ModifyGraph lsize=0.1
	AppendToGraph /C = (0,0,0) ach_1_u5
	appendtograph /r /c=(0,0,65535) ACH_1_pockel5
	Label left "I (pA)"
	Label right "V\\BPockels \\M(V)"
	ModifyGraph width=100,height=125
	AppendToLayout/T Graph3
	ModifyLayout left(Graph3)=225,top(Graph3)=240
	
	display ach_1_p1
	ModifyGraph lsize=0.1
	AppendToGraph /C = (0,0,0) ach_1_u1
	ModifyGraph width=210,height=210
	appendtograph /r /c=(0,0,65535) ACH_1_pockel1
	Label left "I (pA)"
	Label right "V\\BPockels \\M(V)"
	ModifyGraph width=100,height=125
	AppendToLayout/T Graph4
	ModifyLayout left(Graph4)=415,top(Graph4)=240
	
	Display/K=0 AmplitudeData
	ModifyGraph mode=3
	duplicate /o Amplitude_SD amplitude_95ci
	amplitude_95ci = 2*Amplitude_SD
	ErrorBars AmplitudeData Y,wave=(amplitude_95ci,amplitude_95ci)
	wavestats /q amplitudedata
	SetAxis left min(0,wavemin(amplitudedata)),(1.1*v_max)
	SetAxis/A/R left;DelayUpdate
	SetAxis/A/R bottom
	Label left "I (pA)"
	Label bottom "pulse #"
	ModifyGraph width=100,height=125
	AppendToLayout/T Graph5
	ModifyLayout left(Graph5)=605,top(Graph5)=240

	ModifyLayout frame=0

	

tilewindows/o=1/w=(0,0,1000,600)
//	SavePICT/p=output_path /ef=2 /O /win=layout0 /E=-8 as s_filename[0,(strsearch(s_filename,".",0)-1)]+".pdf"
//	SavePICT/O/E=-5/B=72/p=output_path/win=graph0 as (s_filename[0,(strsearch(s_filename,".",0)-1)]+"1-spine.png")
//	SavePICT/O/E=-5/B=72/win=graph1/p=output_path as (s_filename[0,(strsearch(s_filename,".",0)-1)]+"2-trace.full.png")
//	SavePICT/O/E=-5/B=72/win=graph2/p=output_path as (s_filename[0,(strsearch(s_filename,".",0)-1)]+"3-max.response.png")
//	SavePICT/O/E=-5/B=72/win=graph3/p=output_path as (s_filename[0,(strsearch(s_filename,".",0)-1)]+"6-amplitude.png")
//	SavePICT/O/E=-5/B=72/win=graph4/p=output_path as (s_filename[0,(strsearch(s_filename,".",0)-1)]+"4-mid.response.png")
//	SavePICT/O/E=-5/B=72/win=graph5/p=output_path as (s_filename[0,(strsearch(s_filename,".",0)-1)]+"5-min.response.png")
endmacro

function save_results()
string /g wave_data_list
string /g s_filename
string output_base_name = indexedfile(output_path,0,".itx")
	save /b/t/o/p=output_path  wavelist("*_u*",";","")+wavelist("*_p*",";","")+wave_data_list+";average_trace;fit_average_trace" as replacestring(".itx", output_base_name, "_out.itx")
	SavePICT/p=output_path /ef=2 /O /win=layout0 /E=-8 as output_base_name[0,(strsearch(output_base_name,".",0)-1)]+".pdf"
	SavePICT/O/E=-5/B=72/p=output_path/win=layout0 as (output_base_name[0,(strsearch(output_base_name,".",0)-1)]+".png")
	SavePICT/O/E=-5/B=72/p=output_path/win=graph0 as (output_base_name[0,(strsearch(output_base_name,".",0)-1)]+"1-spine.png")
	SavePICT/O/E=-5/B=72/win=graph1/p=output_path as (output_base_name[0,(strsearch(output_base_name,".",0)-1)]+"2-trace.full.png")
	SavePICT/O/E=-5/B=72/win=graph2/p=output_path as (output_base_name[0,(strsearch(output_base_name,".",0)-1)]+"3-max.response.png")
	SavePICT/O/E=-5/B=72/win=graph3/p=output_path as (output_base_name[0,(strsearch(output_base_name,".",0)-1)]+"4-mid-response.png")
	SavePICT/O/E=-5/B=72/win=graph4/p=output_path as (output_base_name[0,(strsearch(output_base_name,".",0)-1)]+"5-min.response.png")
	SavePICT/O/E=-5/B=72/win=graph5/p=output_path as (output_base_name[0,(strsearch(output_base_name,".",0)-1)]+"6-amplitude.png")
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

macro clean_up()
killwindow graph5
killwindow graph4
killwindow graph3
killwindow graph2
killwindow graph1
killwindow graph0
killwindow layout0
killwaves ach_1, ach_3
//killdatafolder /z root:
endmacro


menu "macros"
	"UncagingAnalysis/1"
	"save_results/2"
	"clean_up/3"
end

