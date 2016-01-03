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
	display wave1
	ModifyGraph lsize=0.1
	ModifyGraph width=720,height=225
	AppendToLayout/T Graph1
	ModifyLayout left(Graph1)=20,top(Graph1)=300
	display wave1_p10
//	display wave1_p1
	ModifyGraph lsize=0.1
	AppendToGraph /C = (0,0,0) wave1_u10
//	AppendToGraph /C = (0,0,0) wave1_u1
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

endmacro

macro plot_waves_r()
	duplicate /o AmplitudeData tmp
//	tmp = AmplitudeData[numpnts(amplitudedata)-1-p]
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
//	display wave1_p10
	display wave1_p1
	ModifyGraph lsize=0.1
//	AppendToGraph /C = (0,0,0) wave1_u10
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

endmacro
