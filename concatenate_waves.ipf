#pragma rtGlobals=3		// Use modern global access method and strict wave access.
macro concatenate_waves()
	string wave_list
	variable k= 0
	variable n_waves = 0
	string wave_name
	string root:base_name
	//string base_name 
	wave_list = WaveList("W2015*_1", ";", "" )
	n_waves = itemsinlist(wave_list)
	concatenate /np /o wave_list, wave1
	do
		wave_name = StringFromList(k,wave_list)
		killwaves $wave_name
		k += 1 
	while(k < n_waves)

	wave_list = WaveList("W2015*_2", ";", "" )
	concatenate /np /o wave_list, wave2

	k=0
	do
		wave_name = StringFromList(k,wave_list)
		killwaves $wave_name
		k += 1 
	while(k < n_waves)

	base_name = s_fileName
	print base_name
	UncagingAnalysis("wave1",1,0.04)
endmacro