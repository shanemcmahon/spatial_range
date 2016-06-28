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