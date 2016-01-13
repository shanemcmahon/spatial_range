#pragma rtGlobals=3		// Use modern global access method and strict wave access.


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

newpath current_experiment_dir s_path
print indexedfile(current_experiment_dir,-1,"????")
print itemsinlist(indexedfile(current_experiment_dir,-1,"????"))



macro UncagingAnalysis()
	string/g response_data_file_name, response_number_string, new_response_number_string, new_response_data_file_name
	variable/g response_number, new_response_number
	String traceunc = "ach_1"
	String/G tracepower
	tracepower = "ach_3"
	Variable i=0, Td = 0.004, Tr = 0.002, peak_loc = .007, offsetrange = 0.01, peakrange = 0.002,V_FitMaxIters = 100, v_levelx = 0, t_max_response,  fit_action, fitrange=0.04, this_response
	variable max_power_0=10, min_power_0=0, target_response = 11.5
	Variable/G uncgpnt, betterreturn, V_FitError, cursorA, cursorB, last_power
	variable /g max_power, min_power
	Make/O/D W_coef = NaN					//For holding DiffTwoExp coefficients
	Make/O/D W_sigma = NaN					//For holding DiffTwoExp coefficients

	//---------------------------------------------First check / create all waves

	if ( !WaveExists(root:uncage_power_wave))
		response_data_file_name = s_path+s_filename
		response_data_file_name = response_data_file_name[0,(strsearch(response_data_file_name,".itx",0)-1)]
		response_number_string = response_data_file_name[(strsearch(response_data_file_name,"Line",0)+4),Inf]
		response_number = str2num(response_number_string)
		last_power = read_power_from_prm()
		read_write_prm_0_03(last_power,"C:Documents and Settings:shane:My Documents:sm.power.test.prm","C:Documents and Settings:shane:My Documents:sm.updated.protocol.prm")
	endif
	
	if ( WaveExists(root:uncage_power_wave))
		new_response_number = response_number + 1
		new_response_number_string = PadString(new_response_number_string, (5-floor(log(response_number))), char2num("0") ) + num2str(new_response_number)
		new_response_data_file_name = response_data_file_name[0,(strsearch(response_data_file_name,"Line",0)+3)] + new_response_number_string +".itx"
		LoadWave/T new_response_data_file_name
		new_response_number_string =""
		response_data_file_name = new_response_data_file_name
		response_data_file_name = response_data_file_name[0,(strsearch(response_data_file_name,".itx",0)-1)]
		response_number_string = response_data_file_name[(strsearch(response_data_file_name,"Line",0)+4),Inf]
		response_number = str2num(response_number_string)
	endif

	ach_1 = ach_1 * 100
	
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

		
	i = 0; v_levelx = 0
	findlevel /Q/EDGE=1 /R= (v_levelx,) $tracepower, (0.8*wavemax($tracepower))
	uncgpnt = v_levelx
	findlevel /Q/EDGE=2 /R= (v_levelx,) $tracepower, (0.5*wavemax($tracepower))
	InsertPoints numpnts (uncgpntData), 1, uncgpntData
	uncgpntData[numpnts(uncgpntData)-1] = uncgpnt
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
		// if the exponential fit causes an error, then use a generic estimate for initial parameters
		v_fiterror = 0
		W_coef = {mean($traceunc,(uncgpnt - offsetrange),uncgpnt),-1E-12,.04}
	endif
	W_coef = {(W_coef[1]), uncgpnt, W_coef[2], tr, W_coef[0]}
	Display/N=Checking $traceunc;
	SetAxis bottom (uncgpnt-0.01),(uncgpnt+(2*fitrange))
	do_fit(traceunc, (uncgpnt-fitrange/3), (uncgpnt+fitrange), w_coef)
	AppendToGraph /c=(0,0,0) $("fit_"+traceunc)
	fit_action = PopupChoice ()
	// if (fit_action == 4)	//yes good fit
	// insert_fit_coefs(w_coef, amplitudeData, T0Data, TauDecayData, TauRiseData, OffsetData, TimeDiffUncgpnt_Response, Amplitude_SD, uncgpnt, w_sigma)
	// endif
	if (fit_action == 2)		//Nothing there: save as zero for no response
		w_coef = 0
		w_sigma = NaN
		// insert_fit_coefs(w_coef, amplitudeData, T0Data, TauDecayData, TauRiseData, OffsetData, TimeDiffUncgpnt_Response, Amplitude_SD, uncgpnt, w_sigma)
	endif
	if (fit_action == 3)		//Nothing there: noise! Save NaN
		w_coef = NaN
		w_sigma = NaN
		// insert_fit_coefs(w_coef, amplitudeData, T0Data, TauDecayData, TauRiseData, OffsetData, TimeDiffUncgpnt_Response, Amplitude_SD, uncgpnt, w_sigma)
	endif
	if(fit_action != 1)
		insert_fit_coefs(w_coef, amplitudeData, T0Data, TauDecayData, TauRiseData, OffsetData, TimeDiffUncgpnt_Response, Amplitude_SD, uncgpnt, w_sigma)
	endif
	if(fit_action == 1)
		//No, not good fit. Get fit points from cursors
		DoWindow/K Checking
		Duplicate/O/R=((uncgpnt-(0.5*fitrange)),(uncgpnt+(2*fitrange))) $traceunc, tracecopy
		loess /N=13 /ord=2 /dest=tracecopy srcwave=tracecopy
		Display/N=Checking tracecopy
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

		Display/N=Checking $traceunc;
		SetAxis bottom (uncgpnt-0.01),(uncgpnt+(2*fitrange))
		do_fit(traceunc, cursorA, cursorB, w_coef)
		betterreturn = Refit()
		if (betterreturn ==  1)
			insert_fit_coefs(w_coef, amplitudeData, T0Data, TauDecayData, TauRiseData, OffsetData, TimeDiffUncgpnt_Response, Amplitude_SD, uncgpnt, w_sigma)
		else if (betterreturn == 2)//-------------------------------------If fit not good, saves in waves as NaN
			w_coef = NaN
			w_sigma = NaN
			insert_fit_coefs(w_coef, amplitudeData, T0Data, TauDecayData, TauRiseData, OffsetData, TimeDiffUncgpnt_Response, Amplitude_SD, uncgpnt, w_sigma)
		endif
	endif

	DoWindow/K Checking
	if ( !WaveExists(root:uncage_power_wave))
		Make/O/N=0 uncage_power_wave
		max_power = max_power_0
		min_power = min_power_0
	endif
	last_power = read_power_from_prm()
	this_response = abs(w_coef[0])
	print this_response
	if(this_response > 1.1*target_response)
		max_power = last_power
		//read_write_prm_0_02(((max_power+min_power)/2))
		//"C:Documents and Settings:shane:My Documents:sm.power.test.prm"
		read_write_prm_0_03(((max_power+min_power)/2),"C:Documents and Settings:shane:My Documents:sm.updated.protocol.prm","C:Documents and Settings:shane:My Documents:sm.updated.protocol.prm")
		read_write_prm_0_03(((max_power+min_power)/2),"C:Documents and Settings:shane:My Documents:sm.updated.protocol.prm","C:Documents and Settings:shane:My Documents:sm.power.test.prm")
	endif
	if(this_response < 0.9*target_response)
		min_power = last_power
		//read_write_prm_0_02(((max_power+min_power)/2))
		read_write_prm_0_03(((max_power+min_power)/2),"C:Documents and Settings:shane:My Documents:sm.updated.protocol.prm","C:Documents and Settings:shane:My Documents:sm.updated.protocol.prm")
		read_write_prm_0_03(((max_power+min_power)/2),"C:Documents and Settings:shane:My Documents:sm.updated.protocol.prm","C:Documents and Settings:shane:My Documents:sm.power.test.prm")
	endif
	if(abs(this_response-target_response)<=0.1*target_response)
			read_write_prm_0_03(last_power,"C:Documents and Settings:shane:My Documents:sm.uncage.one.line.prm","C:Documents and Settings:shane:My Documents:sm.updated.protocol.prm")
	endif
	//		KillVariables/A;	KillStrings/A;	KillWaves /A
	killwaves ach_1, ach_3
Endmacro

//macro do_uncage_analysis()
//UncagingAnalysis()
//endmacro

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

function read_write_prm_0_02(new_power)
	variable new_power
	variable file_ref_num
	variable out_file
	variable i_ = 0
	string string_in
	string write_buffer
	variable write_block_size = 2^9
	variable write_buffer_start
	variable write_buffer_stop
	variable power_line_start
	variable power_line_stop
	//	print s_filename[0,(strsearch(s_filename,".",0)-1)]+".prm"
	//	open /R  file_ref_num as (s_path + s_filename[0,(strsearch(s_filename,".",0)-1)]+".prm")
	open /R  file_ref_num as "C:Documents and Settings:shane:My Documents:sm.power.test.prm"
	freadline /T="" file_ref_num, string_in
	close(file_ref_num)
	open out_file as "C:Documents and Settings:shane:My Documents:sm.power.test.prm"
	power_line_start = strsearch(string_in,"DAC6 Train Pulse Pot 1",0)
	power_line_stop = strsearch(string_in,num2char(13),power_line_start)
	string power_line_text = string_in[power_line_start,power_line_stop]
	power_line_text = power_line_text[0,strsearch(power_line_text,"=",0)] + num2str(new_power)+num2char(13)
	do
		write_buffer_start = i_*write_block_size
		write_buffer_stop =(i_+1)*write_block_size-1
		if(write_buffer_stop >= power_line_start)
			break
		endif
		write_buffer = string_in[write_buffer_start,write_buffer_stop]
		fprintf out_file, write_buffer 
		i_ = i_+1
	while(strlen(write_buffer) > 0)
	write_buffer = string_in[write_buffer_start,power_line_start-1]
	fprintf out_file, write_buffer 	
	fprintf out_file, power_line_text	
	string_in = string_in[power_line_stop+1, Inf]
	i_=0
	do
		write_buffer_start = i_*write_block_size
		write_buffer_stop =(i_+1)*write_block_size-1
		write_buffer = string_in[write_buffer_start,write_buffer_stop]
		fprintf out_file, write_buffer 
		i_ = i_+1
	while(strlen(write_buffer) > 0)
	close(out_file)
end

macro do_read_write_prm_0_02()
	variable new_power = 2*pi
	read_write_prm_0_02(new_power)
endmacro

function read_write_prm_0_03(new_power,in_file_name,out_file_name)
	string in_file_name, out_file_name
	variable new_power
	variable file_ref_num
	variable out_file
	variable i_ = 0
	string string_in
	string write_buffer
	variable write_block_size = 2^9
	variable write_buffer_start
	variable write_buffer_stop
	variable power_line_start
	variable power_line_stop
	open /R  file_ref_num as in_file_name
	freadline /T="" file_ref_num, string_in
	close(file_ref_num)
	open out_file as out_file_name
	power_line_start = strsearch(string_in,"DAC6 Train Pulse Pot 1",0)
	power_line_stop = strsearch(string_in,num2char(13),power_line_start)
	string power_line_text = string_in[power_line_start,power_line_stop]
	power_line_text = power_line_text[0,strsearch(power_line_text,"=",0)] + num2str(new_power)+num2char(13)
	do
		write_buffer_start = i_*write_block_size
		write_buffer_stop =(i_+1)*write_block_size-1
		if(write_buffer_stop >= power_line_start)
			break
		endif
		write_buffer = string_in[write_buffer_start,write_buffer_stop]
		fprintf out_file, write_buffer 
		i_ = i_+1
	while(strlen(write_buffer) > 0)
	write_buffer = string_in[write_buffer_start,power_line_start-1]
	fprintf out_file, write_buffer 	
	fprintf out_file, power_line_text	
	string_in = string_in[power_line_stop+1, Inf]
	i_=0
	do
		write_buffer_start = i_*write_block_size
		write_buffer_stop =(i_+1)*write_block_size-1
		write_buffer = string_in[write_buffer_start,write_buffer_stop]
		fprintf out_file, write_buffer 
		i_ = i_+1
	while(strlen(write_buffer) > 0)
	close(out_file)
end

function read_power_from_prm()
	variable file_ref_num, power_line_stop,power_line_start, last_power
	string string_in
	open /R  file_ref_num as "C:Documents and Settings:shane:My Documents:sm.power.test.prm"
	freadline /T="" file_ref_num, string_in
	close(file_ref_num)
	power_line_start = strsearch(string_in,"DAC6 Train Pulse Pot 1",0)
	power_line_stop = strsearch(string_in,num2char(13),power_line_start)
	string power_line_text = string_in[power_line_start,power_line_stop]
	last_power = str2num(power_line_text[strsearch(power_line_text,"=",0)+1,inf])
	return last_power
end