#pragma rtGlobals=3		// Use modern global access method and strict wave access.
macro plot_waves()
	newpath /o output_path s_path
	duplicate /o AmplitudeData tmp
	duplicate /o Amplitude_SD tmp2
	tmp = AmplitudeData[numpnts(amplitudedata)-1-p]
	tmp2 = Amplitude_sd[numpnts(amplitudedata)-1-p]
	amplitudedata = -tmp
	Amplitude_sd = tmp2
	Layout/T /P=Landscape
	NewImage/K=0  $StringFromList(0,WaveList("Trigger*",";",""))
	ModifyGraph width=210,height=210
	ModifyGraph width=210,height={Aspect,1}
	AppendToLayout/T Graph0
	ModifyLayout left(Graph0)=0,top(Graph0)=0
	display ach_1
	Label left "I (pA)"
	ModifyGraph lsize=0.1
	ModifyGraph width=720,height=225
	AppendToLayout/T Graph1
	ModifyLayout left(Graph1)=20,top(Graph1)=300
	display ach_1_p10
	ModifyGraph lsize=0.1
	AppendToGraph /C = (0,0,0) ach_1_u10
	ModifyGraph width=210,height=210
	appendtograph /r /c=(0,0,65535) ACH_1_pockel10
	Label left "I (pA)"
	Label right "V\\BPockels \\M(V)"
	AppendToLayout/T Graph2
	ModifyLayout left(Graph2)=240,top(Graph2)=20
	Display/K=0 AmplitudeData
	ModifyGraph mode=3
	duplicate /o Amplitude_SD amplitude_95ci
	amplitude_95ci = 2*Amplitude_SD
	ErrorBars AmplitudeData Y,wave=(amplitude_95ci,amplitude_95ci)
	wavestats /q amplitudedata
	SetAxis left min(0,wavemin(amplitudedata)),(1.1*v_max)
	Label left "I (pA)"
	Label bottom "pulse #"
	ModifyGraph width=210,height=210
	AppendToLayout/T Graph3
	ModifyLayout left(Graph3)=520,top(Graph3)=20
	ModifyLayout frame=0

	display ach_1_p5
	ModifyGraph lsize=0.1
	AppendToGraph /C = (0,0,0) ach_1_u5
	ModifyGraph width=210,height=210
	appendtograph /r /c=(0,0,65535) ACH_1_pockel5
	Label left "I (pA)"
	Label right "V\\BPockels \\M(V)"

	display ach_1_p1
	ModifyGraph lsize=0.1
	AppendToGraph /C = (0,0,0) ach_1_u1
	ModifyGraph width=210,height=210
	appendtograph /r /c=(0,0,65535) ACH_1_pockel1
	Label left "I (pA)"
	Label right "V\\BPockels \\M(V)"


	SavePICT/p=output_path /ef=2 /O /win=layout0 /E=-8 as s_filename[0,(strsearch(s_filename,".",0)-1)]+".pdf"
	SavePICT/O/E=-5/B=72/p=output_path/win=graph0 as (s_filename[0,(strsearch(s_filename,".",0)-1)]+"1-spine.png")
	SavePICT/O/E=-5/B=72/win=graph1/p=output_path as (s_filename[0,(strsearch(s_filename,".",0)-1)]+"2-trace.full.png")
	SavePICT/O/E=-5/B=72/win=graph2/p=output_path as (s_filename[0,(strsearch(s_filename,".",0)-1)]+"3-max.response.png")
	SavePICT/O/E=-5/B=72/win=graph3/p=output_path as (s_filename[0,(strsearch(s_filename,".",0)-1)]+"6-amplitude.png")
	SavePICT/O/E=-5/B=72/win=graph4/p=output_path as (s_filename[0,(strsearch(s_filename,".",0)-1)]+"4-mid.response.png")
	SavePICT/O/E=-5/B=72/win=graph5/p=output_path as (s_filename[0,(strsearch(s_filename,".",0)-1)]+"5-min.response.png")
endmacro
