#pragma rtGlobals=3		// Use modern global access method and strict wave access.


//*******************************************************************************************************************************

Function PopupChoice ()
	Variable Choose = 0
	Prompt Choose, "Is this fit good?", popup, "Yes: Save fit; No response: Save zero"
	DoPrompt "Goodness of Fit", Choose
  if (V_Flag)
    killvariables v_flag
    abort("User Canceled")								// User canceled
  endif
	return choose
End

Function continue_prompt()
	Variable out = 0
	Prompt out, "Continue?", popup, "Continue; Quit"
	DoPrompt "Continue", out
  if (V_Flag)
    killvariables v_flag
		out = 2								// User canceled
	endif
	return out
End

function do_fit(traceunc, xmin, xmax, w_coef)
	wave w_coef
	string traceunc
	variable xmin, xmax
	variable v_abortcode = 0
	try
		FuncFit/N/Q/H="00000" /NTHR=0 DiffTwoExp2 W_coef  $traceunc(xmin,xmax) /D
	catch
		W_coef = NaN
		Duplicate/O/R=(xmin,xmax) $traceunc, $("fit_"+traceunc)
	endtry

end

macro UncagingAnalysis()
	string/g response_data_file_name, response_number_string, new_response_number_string, new_response_data_file_name
	variable response_number, new_response_number
	String traceunc = "ach_1"
 	string /G protocol_dir = s_path
// string /G protocol_dir = "C:Documents and Settings:shane:My Documents:"
	String tracepower
	tracepower = "ach_3"
	Variable i=0, Td = 0.004, Tr = 0.002, V_FitMaxIters = 100, v_levelx = 0,  fit_action, fitrange=0.04, this_response
	variable max_power_0=10, min_power_0=0, target_response = 0.05
	Variable uncgpnt, V_FitError, last_power
	variable /g max_power, min_power, next_power
	Make/O/D W_coef = NaN					//For holding DiffTwoExp coefficients
  string prm_file_name
  variable temp_var
	//---------------------------------------------First check / create all waves

	if ( !WaveExists(root:uncage_power_wave))
		response_data_file_name = s_path+s_filename
		response_data_file_name = response_data_file_name[0,(strsearch(response_data_file_name,".itx",0)-1)]
		response_number_string = response_data_file_name[(strsearch(response_data_file_name,"Line",0)+4),Inf]
		response_number = str2num(response_number_string)
    prm_file_name = response_data_file_name+".prm"
    last_power = read_power_from_prm(prm_file_name)
    read_write_prm(last_power,response_data_file_name+".prm",protocol_dir+"sm.updated.protocol.prm")
	endif

do //do1
  if ( WaveExists(root:uncage_power_wave))
		new_response_number = response_number + 1
		new_response_number_string = PadString(new_response_number_string, (5-floor(log(response_number))), char2num("0") ) + num2str(new_response_number)
		new_response_data_file_name = response_data_file_name[0,(strsearch(response_data_file_name,"Line",0)+3)] + new_response_number_string +".itx"
		LoadWave/T new_response_data_file_name
		new_response_number_string =""
		response_data_file_name = new_response_data_file_name
		response_data_file_name = response_data_file_name[0,(strsearch(response_data_file_name,".itx",0)-1)]
    prm_file_name = response_data_file_name+".prm"
		response_number_string = response_data_file_name[(strsearch(response_data_file_name,"Line",0)+4),Inf]
		response_number = str2num(response_number_string)
	endif

	i = 0; v_levelx = 0
	findlevel /Q/EDGE=1 /R= (v_levelx,) $tracepower, (0.8*wavemax($tracepower))
	uncgpnt = v_levelx
	findlevel /Q/EDGE=2 /R= (v_levelx,) $tracepower, (0.5*wavemax($tracepower))
  K0 =  mean(ach_1,0,uncgpnt)
  duplicate /O $traceunc temp
  temp = (abs($traceunc - K0))
  wavestats /q temp
  K1 = $traceunc(v_maxloc)-k0
  K2 = td
	CurveFit/N/Q/NTHR=0 exp_XOffset $traceunc(v_maxloc,(uncgpnt + fitrange)) /D
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
		W_coef = {mean($traceunc,0,uncgpnt),-1E-12,.04}
	endif
	W_coef = {(W_coef[1]), uncgpnt, W_coef[2], ((v_maxloc-uncgpnt)*0.33), W_coef[0]}
	Display/N=Checking $traceunc;
	SetAxis bottom (uncgpnt-0.01),(uncgpnt+(2*fitrange))
	do_fit(traceunc, (uncgpnt-fitrange/3), (uncgpnt+fitrange), w_coef)
	AppendToGraph /c=(0,0,0) $("fit_"+traceunc)
	fit_action = PopupChoice ()
	if (fit_action == 2)		//Nothing there: save as zero for no response
		w_coef = 0
	endif
	DoWindow/K Checking
	if ( !WaveExists(root:uncage_power_wave))
		Make/O/N=0 uncage_power_wave
		max_power = max_power_0
		min_power = min_power_0
	endif
  last_power = read_power_from_prm(prm_file_name)
	this_response = abs(w_coef[0])
	print this_response


	// calculate power for next uncaging and update protocol
	// should be moved to dedicated function
	// if(this_response > 1.2*target_response)
	// 	max_power = last_power
	// 	next_power = ((max_power+min_power)/2)
	// 	read_write_prm(next_power,prm_file_name,protocol_dir+"sm.updated.protocol.prm")
	// endif
	// if(this_response < 0.8*target_response)
	// 	min_power = last_power
	// 	next_power = ((max_power+min_power)/2)
  //   read_write_prm(next_power,prm_file_name,protocol_dir+"sm.updated.protocol.prm")
	// endif
	// if(abs(this_response-target_response)<=0.2*target_response)
	// next_power = last_power
	// 		read_write_prm(last_power,protocol_dir+"sm.uncage.one.line.prm",protocol_dir+"sm.updated.protocol.prm")
  //     break
	// endif

	if(next_power_fit(this_response, target_response, last_power))
		read_write_prm(next_power,protocol_dir+"sm.uncage.one.line.prm",protocol_dir+"sm.updated.protocol.prm")
		break
	endif
	print next_power
	read_write_prm(next_power,prm_file_name,protocol_dir+"sm.updated.protocol.prm")

  temp_var = continue_prompt()
  if(temp_var == 2)
  read_write_prm(next_power,protocol_dir+"sm.uncage.one.line.prm",protocol_dir+"sm.updated.protocol.prm")
    break
  endif
	killwaves ach_1, ach_3
while(1)//do1
		KillVariables/A;	KillStrings/A;	KillWaves /z/A
Endmacro


macro fit_power()
DeletePoints (numpnts(ACH_3)-1),1, ACH_3
DeletePoints (numpnts(ACH_4)-1),1, ACH_4

CurveFit/M=2/W=0 poly 3, ACH_4/X=ACH_3/D

duplicate fit_ach_4 ach_4_x
ach_4_x = x
display fit_ach_4 vs ach_4_x
		KillVariables/A;	KillStrings/A;	KillWaves /z/A
endmacro



function next_power_bs(this_response, target_response, last_power)
variable this_response, target_response, last_power
wave ach_4_x, fit_ach_4
variable /g max_power, min_power, next_power
if(this_response > 1.2*target_response)
	max_power = last_power
	next_power = ((max_power+min_power)/2)
	return 0
	// read_write_prm(next_power,prm_file_name,protocol_dir+"sm.updated.protocol.prm")
endif
if(this_response < 0.5 * target_response)
	next_power = 1.41421*interp(last_power, ach_4_x, fit_ach_4)
	next_power = interp(next_power, fit_ach_4,  ach_4_x)
	max_power = max(max_power, 1.1*next_power)

	return 0
	// read_write_prm(next_power,prm_file_name,protocol_dir+"sm.updated.protocol.prm")
endif

if(this_response < 0.8*target_response)
	min_power = last_power
	next_power = ((max_power+min_power)/2)
	return 0
	// read_write_prm(next_power,prm_file_name,protocol_dir+"sm.updated.protocol.prm")
endif
if(abs(this_response-target_response)<=0.2*target_response)
next_power = last_power
return 1
		// read_write_prm(last_power,protocol_dir+"sm.uncage.one.line.prm",protocol_dir+"sm.updated.protocol.prm")
		// break
endif
end

function next_power_fit(this_response, target_response, last_power)
variable this_response, target_response, last_power
wave ach_4_x, fit_ach_4
variable /g max_power, min_power, next_power

if(this_response < 0.5 * target_response)
	next_power = 1.41421*interp(last_power, ach_4_x, fit_ach_4)
	next_power = interp(next_power, fit_ach_4,  ach_4_x)
	max_power = max(max_power, 1.1*next_power)

	return 0
endif

if(abs(this_response-target_response)<=0.2*target_response)
next_power = last_power
return 1
endif
next_power = (((interp(last_power, ach_4_x, fit_ach_4))^2*target_response)/this_response)^0.5
	next_power = interp(next_power, fit_ach_4,  ach_4_x)
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

function read_write_prm(new_power,in_file_name,out_file_name)
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

function read_power_from_prm(prm_in)
  string prm_in
	variable file_ref_num, power_line_stop,power_line_start, last_power
  string string_in
	open /R  file_ref_num as prm_in
	freadline /T="" file_ref_num, string_in
	close(file_ref_num)
	power_line_start = strsearch(string_in,"DAC6 Train Pulse Pot 1",0)
	power_line_stop = strsearch(string_in,num2char(13),power_line_start)
	string power_line_text = string_in[power_line_start,power_line_stop]
	last_power = str2num(power_line_text[strsearch(power_line_text,"=",0)+1,inf])
	return last_power
end

macro uncage_pos_cor(x1,dx,y1,dy)
variable x1,dx,y1,dy
variable/G x2,y2
x2 = x1-dx
y2 = y1-dy
print x2,y2
endmacro
