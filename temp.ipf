
################################################################
################################################################
make marked points panels
regexp for uncaging points in triggersync
Point [0-9]* [x,y,z]=
################################################################
################################################################

redimension/n=(256,256) w_waveaverage

make /o /n=(numpnts(wave0)/3) xpos, ypos
xpos = wave0[p*3]
ypos = wave0[p*3+1]
duplicate /o xpos xpos2
duplicate /o ypos ypos2
xpos2 =128 + (xpos-128)/2
ypos2 =128 + (ypos-128)/2

AppendToGraph/w=graph0 /T ypos vs xpos
ModifyGraph /w=graph0 marker=19,msize=1
ModifyGraph /w=graph0 mode=3
ModifyGraph /w=graph0 noLabel=2,axRGB=(65535,65535,65535)
AppendToGraph/w=graph1 /T ypos vs xpos
ModifyGraph /w=graph1 marker=19,msize=1
ModifyGraph /w=graph1 mode=3
ModifyGraph /w=graph1 noLabel=2,axRGB=(65535,65535,65535)
AutoPositionWindow/M=0/R=graph0 graph1
//AppendToGraph/w=graph2 /T ypos2 vs xpos2
//ModifyGraph /w=graph2 marker=19,msize=1
//ModifyGraph /w=graph2 mode=3
ModifyGraph /w=graph2 noLabel=2,axRGB=(65535,65535,65535)
AutoPositionWindow/M=0/R=graph1 graph2

################################################################
################################################################
analyze spatial range
################################################################
################################################################
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
duplicate /o /r=(0,*)(0,(numpnts(rw_uid)-1)) rw2d_fit_amplitude rw2dNormAmp
columnmins(rw2dNormAmp)
rw2dNormAmp = rw2dNormAmp[p][q]/wcolumnmins[q]
rowmeans(rw2dNormAmp,naremove=1)
duplicate /o wrowmeans rwMeanNormAmp
duplicate /o rw2dNormAmp rw2dNormAmpDiff
rw2dNormAmpdiff = rw2dNormAmp[p][q] -  rwMeanNormAmp[p]
duplicate /o rw2dNormAmpDiff rw2dNormAmpDiffSq
rw2dNormAmpDiffSq = rw2dNormAmpDiff^2
rowmeans(rw2dNormAmpDiffSq,naremove=1)
duplicate /o wrowmeans rwSDNormAmp
rwSDNormAmp = wrowmeans^0.5
duplicate /o rw2dNormAmp rwIsNaN
rwIsNaN = !numtype(rw2dNormAmp[p][q])
rowmeans(rwIsNaN)
wrowmeans = wrowmeans*numpnts(rw_uid)
duplicate /o rwSDNormAmp rwSENormAmp
rwSENormAmp = rwSDNormAmp[p][q]/(wrowmeans[p]^0.5)
display rwMeanNormAmp
setscale /p x,1400,-100,rwMeanNormAmp
ModifyGraph mode=3,marker=19;DelayUpdate
ErrorBars/T=0 rwMeanNormAmp Y,wave=(rwSENormAmp,rwSENormAmp)
K0 = 0;K1 = 1;K2 = 400;
CurveFit/G/H="110"/NTHR=0/TBOX=768/K={0} exp_XOffset  rwMeanNormAmp /D
edit rwMeanNormAmp
appendtotable rwSENormAmp
plotcols(rw2dnormamp)
setscale /p x,1400,-100,"nm",rw2dnormamp
AppendToGraph root:fit_rwMeanNormAmp





pathinfo home; save /b/t wavelist("*",";","") as stringfromlist(9,stringfromlist(0,s_path),":") + stringfromlist(10,stringfromlist(0,s_path),":") + stringfromlist(11,stringfromlist(0,s_path),":")+".itx"
save /b/t wavelist("*",";","") as stringfromlist(9,stringfromlist(0,s_path),":") + stringfromlist(10,stringfromlist(0,s_path),":") + stringfromlist(11,stringfromlist(0,s_path),":")+".itx"
save /b/t/p=home wavelist("*",";","") as stringfromlist(9,stringfromlist(0,s_path),":") + stringfromlist(10,stringfromlist(0,s_path),":") + stringfromlist(11,stringfromlist(0,s_path),":")+".itx"

#pragma rtGlobals=3		// Use modern global access method and strict wave access.
//******************************************************************************
//******************************************************************************
//make figure 1
//******************************************************************************
//******************************************************************************

variable vNumStim
variable vUncageSpacing = 100

newimage /n=UncagePositionZoom48Panel SpineImage48x
make /o /n=(numpnts(wave0)/3) xpos, ypos
xpos = wave0[p*3]
ypos = wave0[p*3+1]
duplicate /o xpos xpos2
duplicate /o ypos ypos2
AppendToGraph/w=UncagePositionZoom48Panel /T ypos vs xpos
ModifyGraph /w=UncagePositionZoom48Panel marker=19,msize=1
ModifyGraph /w=UncagePositionZoom48Panel mode=3
ModifyGraph /w=UncagePositionZoom48Panel noLabel=2,axRGB=(65535,65535,65535)

newimage /n=UncagePositionZoom12Panel SpineImage12x
ModifyGraph /w=UncagePositionZoom12Panel noLabel=2,axRGB=(65535,65535,65535)


ModifyGraph /w=UncagePositionZoom48Panel width=205,height=205
ModifyGraph /w=UncagePositionZoom48Panel margin=20
ModifyGraph /w=UncagePositionZoom12Panel width=205,height=205
ModifyGraph /w=UncagePositionZoom12Panel margin=20



display /n=UncageFullTrace w_uncage_response
wavestats w_uncage_response
SetAxis /w=UncageFullTrace left v_min,(v_max+v_sdev)
duplicate /o w_uncage_time wUncageIndicator
wUncageIndicator = (v_max+v_sdev)
appendtograph /w=UncageFullTrace wUncageIndicator vs w_uncage_time
ModifyGraph /w=UncageFullTrace mode(wUncageIndicator)=3,marker(wUncageIndicator)=2;DelayUpdate
ModifyGraph /w=UncageFullTrace rgb(wUncageIndicator)=(0,0,0)
ModifyGraph /w=UncageFullTrace tick=3,nticks(left)=2,noLabel=2;DelayUpdate
ModifyGraph /w=UncageFullTrace axRGB=(65535,65535,65535),tlblRGB(bottom)=(65535,65535,65535);DelayUpdate
ModifyGraph /w=UncageFullTrace alblRGB(bottom)=(65535,65535,65535)
SetDrawEnv /w=UncageFullTrace ycoord= left;SetDrawEnv dash= 0;DelayUpdate;SetDrawEnv xcoord= bottom
DrawLine /w=UncageFullTrace -0.1,v_min,0.4,v_min
SetDrawEnv /w=UncageFullTrace ycoord= left;SetDrawEnv dash= 0;DelayUpdate;SetDrawEnv xcoord= bottom
DrawLine /w=UncageFullTrace -0.1,v_min,-0.1,v_min+10E-12
ModifyGraph /w=UncageFullTrace width=450,height=102.5
ModifyGraph /w=UncageFullTrace margin(left)=20

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

display /n=SpRange10to14 NormAmpMean
setscale /p x,900,-100, "", NormAmpMean
Label /w=SpRange10to14 left "Normalized Amplitude"
Label /w=SpRange10to14 bottom "Distance (nm)"
ModifyGraph /w=SpRange10to14 mode=3,marker=19;DelayUpdate
ErrorBars/w=SpRange10to14 /T=0 NormAmpMean Y,wave=(NormAmpSE,NormAmpSE)
K0 = 0;K1=64 = 1;K2 = 400;
CurveFit/L=64/n/q/H="110"/NTHR=0/K={0} exp_XOffset  NormAmpMean /D
ModifyGraph /w=SpRange10to14 width=190,height=190
ModifyGraph /w=SpRange10to14 margin=55
ModifyGraph /w=SpRange10to14 margin(top)=0,margin(right)=0
ModifyGraph /w= SpRange10to14 lblMargin=10
SetAxis /w=SpRange10to14 left 0,1





newlayout /n=Fig1
modifylayout /w=Fig1 units=0
modifylayout /w=Fig1 frame = 0
appendtolayout UncagePositionZoom12Panel
modifylayout left(UncagePositionZoom12Panel)=0
modifylayout top(UncagePositionZoom12Panel)=0
appendtolayout UncagePositionZoom48Panel
modifylayout /w=Fig1 frame = 0
modifylayout left(UncagePositionZoom48Panel)=245
modifylayout top(UncagePositionZoom48Panel)=0
appendtolayout UncageFullTrace
modifylayout /w=Fig1 frame = 0
modifylayout left(UncageFullTrace)=0
modifylayout top(UncageFullTrace)=245
appendtolayout UncageResponse0
modifylayout /w=Fig1 frame = 0
modifylayout left(UncageResponse0)=0
modifylayout top(UncageResponse0)=367.5
appendtolayout UncageResponse6
modifylayout /w=Fig1 frame = 0
modifylayout left(UncageResponse6)=166.333
modifylayout top(UncageResponse6)=367.5
appendtolayout UncageResponse10
modifylayout /w=Fig1 frame = 0
modifylayout left(UncageResponse10)=332.666
modifylayout top(UncageResponse10)=367.5
appendtolayout SpRange10to14
modifylayout /w=Fig1 frame = 0
modifylayout left(SpRange10to14)=0
modifylayout top(SpRange10to14)=501


//******************************************************************************
//******************************************************************************
//
//******************************************************************************
//******************************************************************************





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
Save/T/B/p=home  wavelist("*",";","")

edit w_amplitude,w_amplitude_0,w_amplitude_1,w_amplitude_1_se,w_amplitude_se,w_decay_time,w_onset_delay,w_rise_time,w_t0,w_y0
