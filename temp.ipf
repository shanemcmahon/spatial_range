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
