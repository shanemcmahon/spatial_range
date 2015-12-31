#pragma rtGlobals=3		// Use modern global access method and strict wave access.


//*******************************************************************************************************************************

Function PopupChoice ()

	String tracename
	Variable Choose = 4
	Variable/G returnvar
	Prompt traceName,"Trace",popup,TraceNameList("",";",1)
	Prompt Choose, "Is this fit good?", popup, "No: Do a Refit; No response: Save zero; Too noisy, save NaN; Yes: Save fit"
	DoPrompt "Goodness of Fit", Choose
	//why all of the if then logic? It makes the code longer, more difficult to read and maintain. Instead, I would use:
	// return Choose

	return choose
End

function insert_fit_coefs(w_coef, amplitudeData, T0Data, TauDecayData, TauRiseData, OffsetData, TimeDiffUncgpnt_Response, Amplitude_SD, Refwave, commentwav, uncgpnt, w_sigma)
	wave w_coef, amplitudeData, t0data, TauDecayData, TauRiseData, OffsetData, TimeDiffUncgpnt_Response, Amplitude_SD, Refwave, commentwav, w_sigma
	variable uncgpnt

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

function do_fit(traceunc, fitrange, w_coef, tracecopy, fitwave)
	wave w_coef
	string traceunc, tracecopy, fitwave
	variable fitrange
	variable/G uncgpnt
	// Variable/G V_FitError
	variable v_abortcode = 0
	Duplicate/O/R=((uncgpnt-(0.5*fitrange)),(uncgpnt+(2*fitrange))) $traceunc, $tracecopy
	try
		duplicate /O /R=((uncgpnt-fitrange/3),(uncgpnt+fitrange)) $traceunc, fit_wave1
		fit_wave1 = mean($tracecopy)
		FuncFit/N/Q/H="00000" /NTHR=0 DiffTwoExp2 W_coef  $traceunc((uncgpnt-fitrange/3),(uncgpnt+fitrange)) /D
		AppendToGraph /C = (0,0,0) fit_wave1;
		Duplicate/O/R=(uncgpnt,(uncgpnt+fitrange)) fit_wave1, $fitwave
		Duplicate/O fit_wave1, $fitwave
	catch
		W_coef = NaN
		Duplicate/O/R=((uncgpnt-fitrange/3),(uncgpnt+fitrange)) $traceunc, $fitwave
		AppendToGraph /C = (0,0,0) $fitwave;
	endtry

end

Macro UncagingAnalysis (traceunc, uncageinterval, fitrange)
	String traceunc
	Variable uncageinterval=1, fitrange=0.04
	Prompt traceunc, "Response Trace", popup, WaveList("*", ";", "")
	// Prompt uncageinterval, "Interval between uncaging points (s)"
	Prompt fitrange, "Fitting range (s)"

	String/G tracepower, new, newname		//Corresponding power trace saved by PrairieView as _2, while recording as _1
	new = RemoveEnding(traceunc)
	tracepower = new + "2"
	Variable returnchoice
	Variable i=0
	Variable Td = 0.004							//Tau decay estimate
	Variable Tr = 0.002						//Tau rise estimate
	Variable peak_loc = .01
	Variable offsetrange = 0.01					//10ms averaged for offset of trace
	Variable peakrange = 0.001				//1ms averaged for offset of peak
	Variable mini, maxi							//x point range for getting initial valuation of peak position
	Variable offsetpt1, offsetpt2				//Specified later: is the x value range for trace offset
	Variable minpeak, maxpeak					//Peak evaluation range for getting amplitude
	Variable/G amplit						//Variable to store amplitudes into wave
	Variable/G offsetfunction					//w[4] for offsetting trace in DiffTwoExp function
	Variable/G uncgpnt						//uncaging point and it's end
	Variable maximum							//maximum of each uncaging power pulse from trace
	Variable/G betterreturn
	Variable V_FitMaxIters = 100			//Number of iterations FitFunc can go through
	Variable/G Threshold
	Variable/G V_FitError					//If there's a problem e.g. Singular Matrix Error, then saves specifies action (save NaN)
	Variable/G cursorA, cursorB
	Variable W_med
	variable w_mad
	Make/O/D W_coef = NaN					//For holding DiffTwoExp coefficients
	Make/O/D W_sigma = NaN					//For holding DiffTwoExp coefficients
	variable v_levelx = 0
	variable y_range = 0
	variable last_uncage_time = 0
	variable t_max_response
	variable fit_action

	//---------------------------------------------First check / create all waves

	if (!WaveExists(root:Refwave))
		Make/O/T/N=0 Refwave
	endif


	if (!WaveExists(root:amplitudedata))
		Make/O/N=0 AmplitudeData
	endif


	if (!WaveExists(root:T0Data))
		Make/O/N=0 T0Data
	endif


	if (!WaveExists(root:TauDecayData))
		Make/O/N=0 TauDecayData
	endif


	if (!WaveExists(root:TauRiseData))
		Make/O/N=0 TauRiseData
	endif


	if (!WaveExists(root:uncgpntData))
		Make/O/N=0 uncgpntData
	endif


	if (!WaveExists(root:TimeDiffUncgpnt_Response))
		Make/O/N=0 TimeDiffUncgpnt_Response
	endif


	if ( !WaveExists(root:OffsetData))		//this is y0
		Make/O/N=0 OffsetData
	endif


	if ( !WaveExists(root:Amplitude_SD))
		Make/O/N=0 Amplitude_SD
	endif


	if ( !WaveExists(root:commentwav))
		Make/O/T/N=0 commentwav
	endif

	//--------------------------------Start loop
	// Initial loop through the data to find an appropriate y-range for displaying plots
	i = 0
	Duplicate/O $traceunc, smoothed
	Smooth/B=3 1, smoothed
	findlevel /Q /EDGE=1 /R= (v_levelx,) $tracepower, (0.8*wavemax($tracepower))
	findlevel /Q/EDGE=2 /R= (v_levelx,) $tracepower, (0.5*wavemax($tracepower))
	findlevel /Q/EDGE=1 /R= (v_levelx,) $tracepower, (0.8*wavemax($tracepower))
	last_uncage_time = v_levelx
	y_range = wavemax(smoothed,0,v_levelx)-wavemin(smoothed,0,v_levelx)
	do
		findlevel /Q/EDGE=2 /R= (v_levelx,) $tracepower, (0.5*wavemax($tracepower))
		threshold = 0.8*maximum
		maximum = WaveMax($tracepower, v_levelx, +inf)
		if (i >= 1)
			if	(maximum < threshold)
				break
			endif
		endif

		findlevel /Q/EDGE=1 /R= (v_levelx,) $tracepower, (0.8*wavemax($tracepower))
		y_range = max(y_range,wavemax(smoothed,last_uncage_time,v_levelx)-wavemin(smoothed,last_uncage_time,v_levelx))
		last_uncage_time = v_levelx
		i = i + 1
	while(1)

	i = 0; v_levelx = 0
	do
		V_FitError = 0;
		threshold = 0.8*maximum
		maximum = WaveMax($tracepower, v_levelx, +inf)
		if (i >= 1)
			if	(maximum < threshold)
				break
			endif
		endif



		findlevel /Q/EDGE=1 /R= (v_levelx,) $tracepower, (0.8*wavemax($tracepower))
		uncgpnt = v_levelx
		findlevel /Q/EDGE=2 /R= (v_levelx,) $tracepower, (0.5*wavemax($tracepower))
		InsertPoints numpnts (uncgpntData), 1, uncgpntData
		uncgpntData[numpnts(uncgpntData)] = uncgpnt
		K0 = mean($traceunc,(uncgpnt - offsetrange),uncgpnt);		K1 = mean($traceunc,(uncgpnt + peak_loc),(uncgpnt + peak_loc * peakrange)) - k0;K2 = td;
		CurveFit/N/Q/NTHR=0 exp_XOffset $traceunc(uncgpnt + peak_loc,(uncgpnt + peak_loc + fitrange)) /D
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
			v_fiterror = 0
			W_coef = {-1E-12,-1E-12,.04}
		endif
		W_coef = {(W_coef[1]), uncgpnt, W_coef[2], tr, W_coef[0]}; duplicate /o w_coef coef0


		//-------------------------------------Name graphs and Display

		String/G fitwave = nameofwave($traceunc) + "_u" + num2str(i+1)
		String/G tracecopy = nameofwave($traceunc) + "_P" + num2str(i+1)
		Display/N=Checking $traceunc;
		SetAxis/W=Checking left (wavemax($traceunc,(uncgpnt-0.01),(uncgpnt+(2*fitrange)))-y_range), wavemax($traceunc,(uncgpnt-0.01),(uncgpnt+(2*fitrange)));
		SetAxis bottom (uncgpnt-0.01),(uncgpnt+(2*fitrange))
		// do_fit(traceunc, fitrange, w_coef)
		do_fit(traceunc, fitrange, w_coef, tracecopy, fitwave)

		// if ( V_FitError)
		// W_coef = NaN
		// endif

		//---------------------------------------------Start of if-loop popup menu: whether to save data=good fits or NaN=bad fits
		//---------------------------------------------
		fit_action = PopupChoice ()
		if (fit_action == 4)	//yes good fit
			insert_fit_coefs(w_coef, amplitudeData, T0Data, TauDecayData, TauRiseData, OffsetData, TimeDiffUncgpnt_Response, Amplitude_SD, Refwave, commentwav, uncgpnt, w_sigma)
		endif
		if (fit_action == 2)		//Nothing there: save as zero for no response
			w_coef = 0
			w_sigma = NaN
			insert_fit_coefs(w_coef, amplitudeData, T0Data, TauDecayData, TauRiseData, OffsetData, TimeDiffUncgpnt_Response, Amplitude_SD, Refwave, commentwav, uncgpnt, w_sigma)
		endif
		if (fit_action == 3)		//Nothing there: noise! Save NaN
			w_coef = NaN
			w_sigma = NaN
			insert_fit_coefs(w_coef, amplitudeData, T0Data, TauDecayData, TauRiseData, OffsetData, TimeDiffUncgpnt_Response, Amplitude_SD, Refwave, commentwav, uncgpnt, w_sigma)
		endif
		if(fit_action == 1)
			//No, not good fit. Get fit points from cursors
			DoWindow/K Checking
			Duplicate/O/R=((uncgpnt-(0.5*fitrange)),(uncgpnt+(2*fitrange))) $traceunc, $tracecopy
			Duplicate/O/R=(uncgpnt,(uncgpnt+fitrange)) fit_wave1, $fitwave
			loess /N=13 /ord=2 /dest=$tracecopy srcwave=$tracecopy
			Display/N=Checking $tracecopy;	AppendToGraph /C = (0,0,0)$fitwave;	SetAxis/W=Checking left (wavemax($tracecopy)-y_range), (wavemax($tracecopy))
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

			// V_FitError = 0
			// Duplicate/O /R=(cursorA, cursorB) $traceunc, fit
			FuncFit/N/Q/H="00000" /NTHR=0 DiffTwoExp2 W_coef  $traceunc(cursorA, cursorB) /D
			Duplicate/O fit_wave1 fit

			// if ( V_FitError)
			// W_coef = NaN
			// endif
			//---------------------------------------------Display and check if good fit
			String/G fitwave = nameofwave($traceunc) + "_u" + num2str(i+1)
			String/G tracecopy = nameofwave($traceunc) + "_P" + num2str(i+1)
			Duplicate/O/R=(cursorA, cursorB) fit, $fitwave
			Duplicate/O/R=((uncgpnt-(0.5*fitrange)),(uncgpnt+(2*fitrange))) $traceunc, $tracecopy

			Display/N=Checking $tracecopy;	AppendToGraph /C = (0,0,0)$fitwave
			SetAxis/W=Checking left (wavemax($tracecopy)-y_range), (wavemax($tracecopy)); SetAxis/W=Checking bottom (uncgpnt-(0.5*fitrange)),(uncgpnt+(2*fitrange))

			betterreturn = Refit()

			if (betterreturn ==  1)
				insert_fit_coefs(w_coef, amplitudeData, T0Data, TauDecayData, TauRiseData, OffsetData, TimeDiffUncgpnt_Response, Amplitude_SD, Refwave, commentwav, uncgpnt, w_sigma)
			else if (betterreturn == 2)//-------------------------------------If fit not good, saves in waves as NaN
				w_coef = NaN
				w_sigma = NaN
				insert_fit_coefs(w_coef, amplitudeData, T0Data, TauDecayData, TauRiseData, OffsetData, TimeDiffUncgpnt_Response, Amplitude_SD, Refwave, commentwav, uncgpnt, w_sigma)
			endif
		endif

		DoWindow/K Checking
		i += 1
	while(1)
	 Save/T AmplitudeData, T0Data, TauDecayData, TauRiseData, TimeDiffUncgpnt_Response, OffsetData, Amplitude_SD as replacestring(".itx", S_fileName, "_out.itx")
	// KillWaves W_coef, smoothed, fit
	// KillVariables/A;	KillStrings/A;	KillWaves /A
EndMacro

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
