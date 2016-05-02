macro Do_Stim_Response_Average()
	String data_wave_list = s_wavenames
	stim_response_average(data_wave_list)
	//Kill_Wave_List(data_wave_list)
endmacro //Do_Stim_Response_Average

Function stim_response_average(data_wave_list)
//Stim_response_average is a special use function that preprocesses uncaging responses with a particular uncaging pulse design pattern and writes an output file that can be accepted by the uncaging anlysis routine
//The function takes inputs data_wave_list which should contain the names of the uncaging response wave and power wave
//The user specefies the form of the uncaging pulse design pattern with appropriately designed w_refs and w_stimi
//each w_refs entry corresponds to an uncaging event at a specific location, each w_stimi element specifies a pulse number where the specified pulse was deliviered at location i
// for example

// make /wave /o /n=5 w_refs
// make /o/n=5 w_stim1 = {2,3,5,8,12}
// w_refs[0] = w_stim1
// make /o/n=4 w_stim2 = {4,6,9,13}
// w_refs[1] = w_stim2
// make /o/n=3 w_stim3 = {7,10,14}
// w_refs[2] = w_stim3
// make /o/n=2 w_stim4 = {11,15}
// w_refs[3] = w_stim4
// make /o/n=2 w_stim5 = {1,16}
// w_refs[4] = w_stim5

//This specifies at protocol with 5 unique uncaging locations
 // location 1 recieved the 2nd,3rd,5th,8th,and 12th uncaging stimuli
 //location 2 recieved the 4th,6th,9th,and 13th stimuli
 //note stimulus event counting starts at 1, that is the first laser pulse is refered to as stimulus 1, not stimulus 0
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

// loop over uncaging locations
for(j = 0;j<numpnts(w_refs);j+=1)//for2
n_reps = numpnts(w_refs[j])
wave w = w_refs[j]

// for each uncaging location, average the responses
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

end //stim_response_average
