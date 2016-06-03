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