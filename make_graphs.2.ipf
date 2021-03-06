#pragma rtGlobals=3		// Use modern global access method and strict wave access.
macro plot_waves()
	duplicate /o AmplitudeData tmp
	tmp = AmplitudeData[numpnts(amplitudedata)-1-p]
	amplitudedata = -tmp
	Layout/T /P=Landscape
	//	NewImage/K=0  $StringFromList(0,WaveList("Trigger*",";",""))
	//	ModifyGraph width=210,height=210
	ModifyGraph width=210,height={Aspect,1}
	AppendToLayout/T Graph0
	ModifyLayout left(Graph0)=0,top(Graph0)=0
	display ach_1
	ModifyGraph lsize=0.1
	ModifyGraph width=720,height=225
	AppendToLayout/T Graph1
	ModifyLayout left(Graph1)=20,top(Graph1)=300
	//	duplicate /o ach_1_p1 time_
	//	time_ = x
	//	display ach_1_p10 vs time_
	display ach_1_p10
	ModifyGraph lsize=0.1
	AppendToGraph /C = (0,0,0) ach_1_u10
	ModifyGraph width=210,height=210
	SetDrawEnv xcoord= bottom;SetDrawEnv dash= 3;DelayUpdate
	DrawLine t0data[9],0,t0data[9],1
	AppendToLayout/T Graph2
	ModifyLayout left(Graph2)=240,top(Graph2)=20
	Display/K=0 AmplitudeData
	SetAxis left min(0,wavemin(amplitudedata)),*
	ModifyGraph width=210,height=210
	AppendToLayout/T Graph3
	ModifyLayout left(Graph3)=520,top(Graph3)=20
	ModifyLayout frame=0

	display ach_1_p5
	ModifyGraph lsize=0.1
	AppendToGraph /C = (0,0,0) ach_1_u5
	ModifyGraph width=210,height=210
	SetDrawEnv xcoord= bottom;SetDrawEnv dash= 3;DelayUpdate
	DrawLine t0data[4],0,t0data[4],1
	
	display ach_1_p1
	ModifyGraph lsize=0.1
	AppendToGraph /C = (0,0,0) ach_1_u1
	ModifyGraph width=210,height=210
	SetDrawEnv xcoord= bottom;SetDrawEnv dash= 3;DelayUpdate
	DrawLine t0data[0],0,t0data[0],1


	SavePICT /ef=2 /O /win=layout0 /E=-8 as s_filename[0,(strsearch(s_filename,".",0)-1)]+".pdf"
	SavePICT/O/E=-5/B=72/win=graph0 as (s_filename[0,(strsearch(s_filename,".",0)-1)]+"1-spine")
	SavePICT/O/E=-5/B=72/win=graph1 as (s_filename[0,(strsearch(s_filename,".",0)-1)]+"2-trace.full")
	SavePICT/O/E=-5/B=72/win=graph2 as (s_filename[0,(strsearch(s_filename,".",0)-1)]+"3-max.response")
	SavePICT/O/E=-5/B=72/win=graph3 as (s_filename[0,(strsearch(s_filename,".",0)-1)]+"6-amplitude")
	SavePICT/O/E=-5/B=72/win=graph4 as (s_filename[0,(strsearch(s_filename,".",0)-1)]+"4-mid.response")
	SavePICT/O/E=-5/B=72/win=graph5 as (s_filename[0,(strsearch(s_filename,".",0)-1)]+"5-min.response")
endmacro

macro plot_waves_r()
	duplicate /o AmplitudeData tmp
	amplitudedata = -tmp
	Layout/T /P=Landscape
	//	NewImage/K=0  $StringFromList(0,WaveList("Trigger*",";",""))
	//	ModifyGraph width=210,height=210
	ModifyGraph width=210,height={Aspect,1}
	AppendToLayout/T Graph0
	ModifyLayout left(Graph0)=0,top(Graph0)=0
	display wave1
	ModifyGraph lsize=0.1
	ModifyGraph width=720,height=225
	AppendToLayout/T Graph1
	ModifyLayout left(Graph1)=20,top(Graph1)=300
	display wave1_p1
	ModifyGraph lsize=0.1
	AppendToGraph /C = (0,0,0) wave1_u1
	ModifyGraph width=210,height=210
	AppendToLayout/T Graph2
	ModifyLayout left(Graph2)=240,top(Graph2)=20
	Display/K=0 AmplitudeData
	SetAxis left min(0,wavemin(amplitudedata)),*
	ModifyGraph width=210,height=210
	AppendToLayout/T Graph3
	ModifyLayout left(Graph3)=520,top(Graph3)=20
	ModifyLayout frame=0
	SavePICT /ef=2 /O /win=layout0 /E=-8 as s_filename[0,(strsearch(s_filename,".",0)-1)]+".pdf"
	SavePICT/O/E=-5/B=72/win=graph0 as (s_filename[0,(strsearch(s_filename,".",0)-1)]+"1-spine")
	SavePICT/O/E=-5/B=72/win=graph1 as (s_filename[0,(strsearch(s_filename,".",0)-1)]+"2-trace.full")
	SavePICT/O/E=-5/B=72/win=graph2 as (s_filename[0,(strsearch(s_filename,".",0)-1)]+"3-max.response")
	SavePICT/O/E=-5/B=72/win=graph3 as (s_filename[0,(strsearch(s_filename,".",0)-1)]+"4-amplitude")
endmacro
