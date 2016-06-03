#pragma rtGlobals=3		// Use modern global access method and strict wave access.
print s_wavenames
duplicate /o $s_wavenames SpineImage
newimage SpineImage
wave0 = wave0/256*dimsize(SpineImage,0)
wave1 = wave1/256*dimsize(SpineImage,1)
display wave0 vs wave1
AppendToGraph/T wave1 vs wave0
ModifyGraph mode=3


string OutputWaveList = rw_uid;rw2d_response;rw2d_uncage_time;rw2d_fit_start_time;rw2d_fit_stop_time;rw2d_fit_amplitude;rw2d_fit_amplitude_se;rw2d_fit_t0;rw2d_fit_decay_time;rw2d_fit_rise_time;rw2d_fit_y0;
rw2d_fit_onset_delay;rw2d_amplitude_0;rw2d_amplitude_0_np;rw2dFitAmplitude1;rw2dFitAmplitudeSE0;rw2dFitAmplitudeSE1;rw2dFitAmplitude0;rw3d_uncaging_response;rw3d_fits;W2dNrAmplitude0;W2dNrAmplitude1;
W2dNrAmplitude2;WAmplitudeCorrelation;WAmplitudeNrCorrelation;rwPockelsVoltage;
save /b /o OutputWaveList

print (wfitamplitude1[2]-wfitamplitude[2])/(wfitamplitude1SE[2]^2+wfitamplitudese[2]^2)^0.5


duplicate /o /r= [0,9][2] w2dFitUncagingResponseCoef, wtemp
wtemp = wtemp - wAvgUncageResponseFitCoef[2]
duplicate /o wtemp wtemp2
wtemp2 = W2dfitUncagingResponseCoefSE[p][2]^2 + WAvgUncageResponseFitCoefSE[2]^2
wtemp2 = wtemp2^0.5
duplicate /o wtemp2 wtemp3
wtemp3 = wtemp/wtemp2


duplicate /o /r= [0,9][3] w2dFitUncagingResponseCoef, wtemp
wtemp = wtemp - wAvgUncageResponseFitCoef[3]
duplicate /o wtemp wtemp2
wtemp2 = W2dfitUncagingResponseCoefSE[p][3]^2 + WAvgUncageResponseFitCoefSE[3]^2
wtemp2 = wtemp2^0.5
duplicate /o wtemp2 wtemp4
wtemp4 = wtemp/wtemp2

duplicate /o w_amplitude_1 wAmpDiff
wAmpDiff=w_amplitude_1 - w_amplitude
duplicate /o w_amplitude_1 wAmpDiffSE
wAmpDiffSE = (w_amplitude_1_se^2 + w_amplitude_se^2)^0.5
duplicate/o w_amplitude_1 wAmpDiffZ
wAmpDiffZ = wAmpDiff/wAmpDiffSE

//edit wtemp3,wtemp4,wAmpDiffZ,w_amplitude_0,w_amplitude_1,w_amplitude


duplicate /o w_avg_response wAvgResponseResid
wAvgResponseResid = w_avg_response - fit_w_avg_response

duplicate /o w_avg_response wResponseTime
wResponseTime = x
duplicate /o wResponseTime fit_w_avg_response
fit_w_avg_response = difftwoexp2(wAvgUncageResponseFitCoef,wResponseTime)
wAvgResponseResid = w_avg_response - fit_w_avg_response
statsresample /N =(numpnts(w_avg_response)) wAvgResponseResid
