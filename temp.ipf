#pragma rtGlobals=3		// Use modern global access method and strict wave access.
variable vNumStim
variable vUncageSpacing = 100

dowindow/k graph0
display w_uncage_response
wavestats w_uncage_response
SetAxis left v_min,(v_max+v_sdev)
duplicate /o w_uncage_time wUncageIndicator
wUncageIndicator = (v_max+v_sdev)
appendtograph wUncageIndicator vs w_uncage_time
ModifyGraph mode(wUncageIndicator)=3,marker(wUncageIndicator)=2;DelayUpdate
ModifyGraph rgb(wUncageIndicator)=(0,0,0)
ModifyGraph width={Aspect,7}
ModifyGraph tick=3,nticks(left)=2,noLabel=2;DelayUpdate
ModifyGraph manTick(left)={0,10,-12,0},manMinor(left)={0,50};DelayUpdate
ModifyGraph manTick(bottom)={2,1,0,0},manMinor(bottom)={0,0};DelayUpdate
ModifyGraph axRGB=(65535,65535,65535),tlblRGB(bottom)=(65535,65535,65535);DelayUpdate
ModifyGraph alblRGB(bottom)=(65535,65535,65535)
SetDrawEnv ycoord= left;SetDrawEnv dash= 0;DelayUpdate;SetDrawEnv xcoord= bottom
DrawLine -0.1,v_min,0.4,v_min
SetDrawEnv ycoord= left;SetDrawEnv dash= 0;DelayUpdate;SetDrawEnv xcoord= bottom
DrawLine -0.1,v_min,-0.1,v_min+10E-12
ModifyGraph width=498.898,height=113.386
ModifyGraph margin=10

vNumStim = dimsize(w2d_responses,0)
duplicate /o w_uncage_time wUncageIndicatorTime
wUncageIndicatorTime = w_uncage_time[0]-w_fit_start_time[0]

dowindow/k graph1
Display w2d_responses[0][*]
appendtograph w2d_fits[0][*]
ModifyGraph rgb(w2d_fits)=(0,0,0)
appendtograph wUncageIndicator vs wUncageIndicatorTime
ModifyGraph mode(wUncageIndicator)=3,marker(wUncageIndicator)=2;DelayUpdate
ModifyGraph rgb(wUncageIndicator)=(0,0,0)
SetAxis left v_min,(v_max+v_sdev)
ModifyGraph tick=3,nticks(left)=2,noLabel=2;DelayUpdate
ModifyGraph manTick(left)={0,10,-12,0},manMinor(left)={0,50};DelayUpdate
ModifyGraph manTick(bottom)={2,1,0,0},manMinor(bottom)={0,0};DelayUpdate
ModifyGraph axRGB=(65535,65535,65535),tlblRGB(bottom)=(65535,65535,65535);DelayUpdate
ModifyGraph alblRGB(bottom)=(65535,65535,65535)
SetDrawEnv ycoord= left;SetDrawEnv dash= 0;DelayUpdate;SetDrawEnv xcoord= bottom
DrawLine -0.005,v_min,0.0,v_min
SetDrawEnv ycoord= left;SetDrawEnv dash= 0;DelayUpdate;SetDrawEnv xcoord= bottom
DrawLine -0.005,v_min,-0.005,v_min+10E-12
ModifyGraph margin=10,width=136.2993,height=113.386
ModifyGraph margin(left)=20


dowindow/k graph2
Display w2d_responses[5][*]
appendtograph w2d_fits[5][*]
ModifyGraph rgb(w2d_fits)=(0,0,0)
appendtograph wUncageIndicator vs wUncageIndicatorTime
ModifyGraph mode(wUncageIndicator)=3,marker(wUncageIndicator)=2;DelayUpdate
ModifyGraph rgb(wUncageIndicator)=(0,0,0)
SetAxis left v_min,(v_max+v_sdev)
ModifyGraph tick=3,nticks(left)=2,noLabel=2;DelayUpdate
ModifyGraph manTick(left)={0,10,-12,0},manMinor(left)={0,50};DelayUpdate
ModifyGraph manTick(bottom)={2,1,0,0},manMinor(bottom)={0,0};DelayUpdate
ModifyGraph axRGB=(65535,65535,65535),tlblRGB(bottom)=(65535,65535,65535);DelayUpdate
ModifyGraph alblRGB(bottom)=(65535,65535,65535)
SetDrawEnv ycoord= left;SetDrawEnv dash= 0;DelayUpdate;SetDrawEnv xcoord= bottom
DrawLine -0.005,v_min,0.0,v_min
SetDrawEnv ycoord= left;SetDrawEnv dash= 0;DelayUpdate;SetDrawEnv xcoord= bottom
DrawLine -0.005,v_min,-0.005,v_min+10E-12
ModifyGraph margin=10,width=136.2993,height=113.386
ModifyGraph margin(left)=20

dowindow/k graph3
Display w2d_responses[9][*]
appendtograph w2d_fits[9][*]
ModifyGraph rgb(w2d_fits)=(0,0,0)
appendtograph wUncageIndicator vs wUncageIndicatorTime
ModifyGraph mode(wUncageIndicator)=3,marker(wUncageIndicator)=2;DelayUpdate
ModifyGraph rgb(wUncageIndicator)=(0,0,0)
SetAxis left v_min,(v_max+v_sdev)
ModifyGraph tick=3,nticks(left)=2,noLabel=2;DelayUpdate
ModifyGraph manTick(left)={0,10,-12,0},manMinor(left)={0,50};DelayUpdate
ModifyGraph manTick(bottom)={2,1,0,0},manMinor(bottom)={0,0};DelayUpdate
ModifyGraph axRGB=(65535,65535,65535),tlblRGB(bottom)=(65535,65535,65535);DelayUpdate
ModifyGraph alblRGB(bottom)=(65535,65535,65535)
SetDrawEnv ycoord= left;SetDrawEnv dash= 0;DelayUpdate;SetDrawEnv xcoord= bottom
DrawLine -0.005,v_min,0.0,v_min
SetDrawEnv ycoord= left;SetDrawEnv dash= 0;DelayUpdate;SetDrawEnv xcoord= bottom
DrawLine -0.005,v_min,-0.005,v_min+10E-12
ModifyGraph margin=10,width=136.2993,height=113.386
ModifyGraph margin(left)=20



dowindow/k graph4
display w_amplitude
setscale /p x,((numpnts(w_amplitude)-1)*vUncageSpacing),-vUncageSpacing,"nm",w_amplitude
Label bottom "distance \\u#2 (nm)"
ModifyGraph mode=3,marker=19;DelayUpdate
ErrorBars/T=0/L=0.7 w_amplitude Y,wave=(w_amplitude_se,w_amplitude_se)
Label left "I\\u#2 (pA)"
SetAxis left *,max(wavemax(w_amplitude),0)









duplicate /o  rw2d_fit_amplitude rw2dFitAmplitude
Redimension /n=(-1,numpnts(rw_uid)) rw2dFitAmplitude
columnmins(rw2dFitAmplitude)
duplicate /o wcolumnmins wMaxAmp
duplicate /o rw2dFitAmplitude rw2dNormAmplitude
rw2dNormAmplitude = rw2dFitAmplitude[p][q]/wMaxAmp[q]
rowmeans(rw2dNormAmplitude,naremove=1)

duplicate /o  rw2d_response rw2dResponse
Redimension /n=(-1,numpnts(rw_uid)) rw2dResponse
rowmeans(rw2dResponse)
rename wrowmeans ResponseTraceAverage

duplicate /o  rw2d_fit_amplitude rw2dFitAmplitude
Redimension /n=(-1,numpnts(rw_uid)) rw2dFitAmplitude
rowmeans(rw2dfitamplitude,naremove=1)
rename wrowmeans wFitAmplitudeMean


variable Rs
variable i0
variable i1
variable dv
variable Ri
variable Vr
string RsTraceName
RsTraceName = "root:W20160614g5_5_12_1_1"
duplicate /o $RsTraceName RsTrace
duplicate /o $(RsTraceName+"_Amp") RsTraceAmp
duplicate /o $(RsTraceName+"_Dur") RsTraceDur

string RiTraceName
RiTraceName = "root:W20160614g5_5_13_1_1"
duplicate /o $RiTraceName RiTrace
duplicate /o $(RiTraceName+"_Amp") RiTraceAmp
duplicate /o $(RiTraceName+"_Dur") RiTraceDur

Rs = (RsTraceAmp[2] - RsTraceAmp[3])/(mean(RsTrace,RsTraceDur[1],RsTraceDur[2])-wavemin(RsTrace))
print Rs
i0 = mean(RiTrace,RiTraceDur[1],RiTraceDur[2])
i1 = mean(RiTrace,RiTraceDur[4]-.25,RiTraceDur[4])
dv = RiTraceAmp[3] - RiTraceAmp[2]
Ri = dv/(i1-i0)
Vr = RiTraceAmp[0]-Ri*i0
print Vr
print Rs,Ri,Vr


Save/T/B  wavelist("ws*1p1",";","")+wavelist("ws*2p1",";","") as wavelist("ws*1p1","","")+"++.itx"

edit w_amplitude,w_amplitude_0,w_amplitude_1,w_amplitude_1_se,w_amplitude_se,w_decay_time,w_onset_delay,w_rise_time,w_t0,w_y0
