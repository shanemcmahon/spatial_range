#pragma rtGlobals=3
macro do_analyze_response()
string /g s_path
if(!exists("ach_1"))
if(exists("s_path"))
newpath /o wd, s_path
LoadWave/T /p=wd
else
newpath /o wd
LoadWave/T /p=wd
endif
endif
analyze_response()
endmacro

function analyze_response()
variable uncgpnt, peak_time
variable td = 0.005
wave ach_1, ach_3
variable target_response = 8
uncgpnt = 9.5

prompt target_response, "Target response amplitude (pA)"
prompt uncgpnt, "Uncage pulse start time"
doprompt "Set variables",target_response, uncgpnt

duplicate /o /r=(uncgpnt-0.05,uncgpnt+0.1) ach_1 last_response_wave
duplicate /o /r=(uncgpnt-0.05,uncgpnt+0.1) ach_3 last_stim_wave

K0 =  mean(last_response_wave,0,uncgpnt)
Loess/pass=1 /ord=1 /n=(2^5-1)/DEST=temp srcWave=last_response_wave
duplicate/o temp temp2
temp2 = (abs(temp2 - K0))
wavestats /q temp2
peak_time = v_maxloc
K1 =  temp(v_maxloc)-k0
K2 = td
//K0 = 0;K1 = 15;K2 = 0.005;
killwaves /z w_sigma
duplicate/o /r=(peak_time,) last_response_wave decay_phase_wave
CurveFit/q/G/NTHR=0/K={peak_time} exp_XOffset  decay_phase_wave /D
make/o /n=5 w_coef
w_coef = {k1,uncgpnt,k2,((peak_time-uncgpnt)/3),k0}
FuncFit/N/Q/H="00000" /NTHR=0 DiffTwoExp2 W_coef  last_response_wave((uncgpnt-k2),(peak_time+3*k2)) /D
Display/K=0 last_response_wave
AppendToGraph  root:fit_last_response_wave
ModifyGraph rgb(fit_last_response_wave)=(0,0,0)

if(!exists("amplitude_data"))
make /n=1 amplitude_data
make /n=1 stim_data
amplitude_data = abs(w_coef[0])
stim_data = mean(last_stim_wave, (uncgpnt+.0001),(uncgpnt+.0005) )
killwaves /z ach_1, ach_3
return 1
endif

InsertPoints 0, 1, amplitude_data, stim_data
amplitude_data[0] = abs(w_coef[0])
stim_data[0] = mean(last_stim_wave, (uncgpnt+.0001),(uncgpnt+.0005) )
killwaves /z ach_1, ach_3
end


menu "macros"
"analyze_response"
"update_prm/2"
"make_power_lut"
"uncage_pos_cor"
"do_analyze_response/1"
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

function next_power_fit(this_response, target_response, last_power)
	variable this_response, target_response, last_power
	variable /g max_power, min_power, next_power
	wave uncage_power_lut, uncage_volt_lut

		// duplicate /o uncage_power_lut temp
		// temp = temp -	uncage_power_lut(last_power)
		// temp = abs(temp)
		// wavestats /q temp
		// last_power = uncage_volt_lut[x2pnt(temp,v_minloc)]

	if(abs(this_response-target_response)<=0.2*target_response)
		next_power = last_power
		return 1
	endif

	if(this_response == 0)
		duplicate /o uncage_power_lut temp
		temp = temp -	uncage_power_lut(last_power)*1.41421
		temp = abs(temp)
		wavestats /q temp
		next_power = uncage_volt_lut[x2pnt(temp, v_minloc)]
		min_power = last_power
		duplicate /o uncage_power_lut temp
		temp = temp -	uncage_power_lut(last_power)*1.41421
		temp = abs(temp)
		wavestats /q temp
		max_power = max(uncage_volt_lut[x2pnt(temp, v_minloc)],max_power)

		return 0
	endif

	next_power = uncage_power_lut(last_power) * ((target_response/this_response)^0.5)
	duplicate /o uncage_power_lut temp
	temp = temp -	next_power
	temp = abs(temp)
	wavestats /q temp
	next_power = uncage_volt_lut[x2pnt(temp, v_minloc)]
	return 0
end

function make_power_lut()
	string wave_list = wavelist("ach_4*",";","")
	variable i, npts
	wave ach_4_cycle1,ach_4_cycle2,ach_4_cycle3,ach_4_cycle4,ach_4_cycle5,ach_4_cycle6,ach_4_cycle7
	npts = numpnts($stringfromlist(0,wave_list))

	ach_4_cycle1[0,1] = ach_4_cycle1[2]
	ach_4_cycle1[npts-1] = ach_4_cycle1[npts-2]
	ach_4_cycle2[0,1] = ach_4_cycle2[2]
	ach_4_cycle2[npts-1] = ach_4_cycle2[npts-2]
	ach_4_cycle3[0,1] = ach_4_cycle3[2]
	ach_4_cycle3[npts-1] = ach_4_cycle3[npts-2]
	ach_4_cycle4[0,1] = ach_4_cycle4[2]
	ach_4_cycle4[npts-1] = ach_4_cycle4[npts-2]
	ach_4_cycle5[0,1] = ach_4_cycle5[2]
	ach_4_cycle5[npts-1] = ach_4_cycle5[npts-2]
	ach_4_cycle6[0,1] = ach_4_cycle6[2]
	ach_4_cycle6[npts-1] = ach_4_cycle6[npts-2]
	ach_4_cycle7[0,1] = ach_4_cycle7[2]
	ach_4_cycle7[npts-1] = ach_4_cycle7[npts-2]
	concatenate /np /o wave_list, uncage_power
	make /n = 140 /o uncage_power_lut
	setscale /p x, 0, 0.05, uncage_power_lut
	for(i = 0;i < 140;i = i+1)	// Initialize variables;continue test
		uncage_power_lut[i] = mean(uncage_power,pnt2x(uncage_power,i*200),pnt2x(uncage_power,(i+1)*200))
	endfor												// Execute body code until continue test is FALSE
	duplicate /o uncage_power_lut uncage_volt_lut
	uncage_volt_lut = x
	display uncage_power_lut vs uncage_volt_lut
	make /n=2 power_0
	Edit/K=0 power_0;DelayUpdate
	killwaves /a/z
end

function update_prm()
variable /g target_response = 8
variable /g max_power, min_power, next_power, last_power
string prm_file_name, protocol_dir, response_data_file_name
variable temp_var, converge_indicator, this_response
svar s_path,s_filename
wave amplitude_data, uncage_power_lut, uncage_volt_lut, stim_data
string path_list =   ("C:Documents and Settings:shane:My Documents:;"+s_path)
temp_var = do_variable_prompt(path_list)
protocol_dir = stringfromlist((temp_var-1), path_list)
response_data_file_name = s_path+s_filename
response_data_file_name = response_data_file_name[0,(strsearch(response_data_file_name,".itx",0)-1)]
prm_file_name = response_data_file_name+".prm"
print prm_file_name

last_power = mean(stim_data)
//duplicate /o uncage_power_lut temp
duplicate /o uncage_volt_lut temp
temp = temp -	last_power
temp = abs(temp)
wavestats /q temp
last_power = uncage_volt_lut[x2pnt(temp,v_minloc)]

read_write_prm(last_power,response_data_file_name+".prm",protocol_dir+"sm.updated.protocol.prm")

this_response = mean(amplitude_data)
print this_response,target_response,last_power
converge_indicator = next_power_fit(this_response, target_response, last_power)
print "Last power",last_power,uncage_power_lut(last_power)
print "Next power",next_power,uncage_power_lut(next_power)

duplicate /o uncage_power_lut temp
temp = temp -	uncage_power_lut(next_power)/1.41421
temp = abs(temp)
wavestats /q temp
read_write_prm(uncage_volt_lut[x2pnt(temp, v_minloc)]	,protocol_dir+"sm.uncage.one.line.prm",protocol_dir+"sm.uncage.line.0.prm")
read_write_prm(next_power,protocol_dir+"sm.uncage.one.line.prm",protocol_dir+"sm.uncage.line.1.prm")
duplicate /o uncage_power_lut temp
temp = temp -	uncage_power_lut(next_power)*1.41421
temp = abs(temp)
wavestats /q temp
read_write_prm(uncage_volt_lut[x2pnt(temp, v_minloc)]	,protocol_dir+"sm.uncage.one.line.prm",protocol_dir+"sm.uncage.line.2.prm")

end

function do_variable_prompt(string_list)
	string string_list
	variable /g target_response
	variable /g max_power_0, min_power_0
	variable  return_var, temp_var, tmp_max_power, tmp_min_power
	wave power_0
	temp_var = target_response
	tmp_max_power = power_0[1]
	tmp_min_power = power_0[0]
	prompt tmp_max_power, "Max power"
	prompt tmp_min_power, "Min power"
	prompt temp_var, "Target response amplitude (pA)"
	Prompt return_var, "Protocol directory", popup, string_list
	doprompt "Set variables",return_var, temp_var, tmp_max_power, tmp_min_power
	target_response = temp_var
	max_power_0 = tmp_max_power
	min_power_0 = tmp_min_power
	return(return_var)
end

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

macro uncage_pos_cor(x1,dx,y1,dy)
	variable x1,dx,y1,dy
	variable/G x2,y2
	x2 = x1-dx
	y2 = y1-dy
	print x2,y2
endmacro

