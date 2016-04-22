#pragma rtGlobals=3		// Use modern global access method and strict wave access.

Function Kill_Wave_List(wave_list)
//
// Kill_Wave_List kills the waves listed in wave_list
//
String wave_list //A string containing a list of wave names
Variable n_items = ItemsInList(wave_list) //variable to hold the number of items in the wave list
Variable i // variable for iteration control
	for(i=0;i<n_items;i+=1)	// for1
		print StringFromList(i, wave_list)
		wave w = $StringFromList(i, wave_list)
		killwaves /z w
	endfor												//for1
end
