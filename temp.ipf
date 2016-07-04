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
Rs = (W20160607g3_3_13_1_1_Amp[2] - W20160607g3_3_13_1_1_Amp[3])/(mean(W20160607g3_3_13_1_1,W20160607g3_3_13_1_1_Dur[1],W20160607g3_3_13_1_1_Dur[2])-wavemin(W20160607g3_3_13_1_1))
print
i0 = mean(W20160607g3_3_14_1_1,W20160607g3_3_14_1_1_Dur[1],W20160607g3_3_14_1_1_Dur[2])
i1 = mean(W20160607g3_3_14_1_1,W20160607g3_3_14_1_1_Dur[4]-.25,W20160607g3_3_14_1_1_Dur[4])
dv = W20160607g3_3_14_1_1_Amp[3] - W20160607g3_3_14_1_1_Amp[2]
Ri = dv/(i1-i0)
