#pragma rtGlobals=3		// Use modern global access method and strict wave access.
//variable tmpv1 = 0
//duplicate /O amplitudeData A
//duplicate /O Amplitude_SD A_SE
//tmpv1 = abs(wavemin(A))
//A = A/tmpv1
tmpv1 = abs(wavemin(AmplitudeData))
AmplitudeData = AmplitudeData/tmpv1
//duplicate /O Amplitudedata tmp; AmplitudeData = tmp[9-p]
concatenate /O {A, AmplitudeData}, tmp
killwaves A; duplicate /O tmp A; killwaves tmp
killwaves amplitudedata,t0data,taudecaydata,taurisedata,timediffuncgpnt_response,offsetdata,amplitude_sd

