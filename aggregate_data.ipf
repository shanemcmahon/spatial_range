function make_folder()
string record_id_l
svar record_id
string path_l
svar path
prompt record_id_l,"record_id"
prompt path_l,"path"
record_id_l=record_id
path_l=path
doprompt "enter variables",record_id_l,path_l
record_id = record_id_l
path = path_l
print record_id, path
setdatafolder $path
newdatafolder/o/s $record_id
end

macro copy_data()
string/g new_name
string/g new_path
new_name = "a_"+root:record_id
duplicate amplitudedata, $new_name
new_path = root:path+":amplitude:"
movewave $new_name, $new_path

new_name = "tr_"+root:record_id
duplicate taurisedata, $new_name
new_path = root:path+":tr:"
movewave $new_name, $new_path


new_name = "td_"+root:record_id
duplicate taudecaydata, $new_name
new_path = root:path+":td:"
movewave $new_name, $new_path


new_name = "dt_"+root:record_id
duplicate timediffuncgpnt_response, $new_name
new_path = root:path+":dt:"
movewave $new_name, $new_path

setdatafolder root:
endmacro

macro align_amplitude_waves()
string waves
variable n_waves,i
waves = wavelist("*",";","")
n_waves =  itemsinlist(waves)
i = 0
do
reverse $StringFromList(i,waves)
InsertPoints 0, (10-numpnts($StringFromList(i,waves))), $StringFromList(i,waves)
i = i+1
while(i<n_waves)
end macro

function average_waves()
string expression = ""
string wave_out
string waves
variable n_waves,i, npts
prompt wave_out,"output wave name"
prompt expression,"wave name filter"
doprompt "",expression
waves = wavelist((expression+"*"),";","")
n_waves =  itemsinlist(waves)
i = 0
npts = (numpnts($StringFromList(0,waves)))
make /o/n=(numpnts($StringFromList(0,waves))) temp_avg
make /o/n=(numpnts($StringFromList(0,waves))) temp_se
killwaves /z wave_2d
concatenate waves, wave_2d
make/o /n=(n_waves) temp
do
temp = wave_2d[i][p]
wavestats /q temp
temp_avg[i] = v_avg
temp_se[i] = v_sem
i = i+1
while(i<npts)
//make /o/n=(numpnts($StringFromList(0,waves))) $"average_"+expression
duplicate /o temp_avg $"average_"+expression
//make /o/n=(numpnts($StringFromList(0,waves))) $"se_"+expression
duplicate /o temp_se $"se_"+expression
end

function normalize_waves()
string expression = ""
string mode
string waves
variable n_waves,i, npts
prompt expression,"wave name filter"
prompt mode,"normalization_mode",popup,"max;mean"
doprompt "",expression,mode
if(stringmatch(mode,"mean"))
waves = wavelist((expression+"*"),";","")
n_waves =  itemsinlist(waves)
i = 0
do
duplicate/o $StringFromList(i,waves) temp
temp = abs(temp)
duplicate/o $StringFromList(i,waves) temp2
temp2 = temp2/mean(temp)
duplicate temp2 $("norm_"+StringFromList(i,waves))
i = i+1
while(i<n_waves)
return 0
endif
if(stringmatch(mode,"max"))
waves = wavelist((expression+"*"),";","")
n_waves =  itemsinlist(waves)
i = 0
do
duplicate/o $StringFromList(i,waves) temp
temp = abs(temp)
duplicate/o $StringFromList(i,waves) temp2
temp2 = temp2/wavemax(temp)
duplicate temp2 $("norm_"+StringFromList(i,waves))
i = i+1
while(i<n_waves)
return 0
endif
abort("1")

end

menu "macros"
"make_folder/4"
"copy_data/5"
"align_amplitude_waves"
"average_waves"
"normalize_waves"
end

K0 = 0;K1 = -1;K2 = 300;
CurveFit/G/H="110"/NTHR=0/K={0} exp_XOffset  average_norm_a_r /D
