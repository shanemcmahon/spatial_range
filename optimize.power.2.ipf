#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
variable uncgpnt, peak_time
variable td = 0.005

duplicate /o /r=(9.45,9.6) ach_1 last_response_wave
uncgpnt = 9.5
K0 =  mean(last_response_wave,0,uncgpnt)
Loess/pass=1 /ord=1 /n=(2^5-1)/DEST=temp srcWave=last_response_wave
duplicate/o temp temp2
temp2 = (abs(temp2 - K0))
wavestats /q temp2
peak_time = v_maxloc
K1 =  temp(v_maxloc)-k0
K2 = td
//K0 = 0;K1 = 15;K2 = 0.005;
killwaves /z w_sigma
duplicate /r=(peak_time,) last_response_wave decay_phase_wave
CurveFit/G/NTHR=0/K={peak_time} exp_XOffset  decay_phase_wave /D
make/o /n=5 w_coef
w_coef = {k1,uncgpnt,k2,((peak_time-uncgpnt)/3),k0}
FuncFit/N/Q/H="00000" /NTHR=0 DiffTwoExp2 W_coef  last_response_wave((uncgpnt-k2),(peak_time+3*k2)) /D
Display/K=0 last_response_wave
AppendToGraph  root:fit_last_response_wave
ModifyGraph rgb(fit_last_response_wave)=(0,0,0)


Function DiffTwoExp2(w,t) : FitFunc
	Wave w
	Variable t

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/    IF    ((t-t0) < 0 )
	//CurveFitDialog/        f(t) = 0
	//CurveFitDialog/    ELSE
	//CurveFitDialog/        f(t) = (gsyn*(-exp(-(t-t0)/tr)+exp(-(t-t0)/td))/(-exp(-(t0+(td*tr/(td-tr))*ln(td/tr)-t0)/tr)+exp(-(t0+(td*tr/(td-tr))*ln(td/tr)-t0)/td)) -y0)
	//CurveFitDialog/    ENDIF
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ t
	//CurveFitDialog/ Coefficients 5
	//CurveFitDialog/ w[0] = gsyn
	//CurveFitDialog/ w[1] = t0
	//CurveFitDialog/ w[2] = td
	//CurveFitDialog/ w[3] = tr
	//CurveFitDialog/ w[4] = y0
	IF    ((t-w[1]) < 0 )
		return w[4]
	ELSE
		return (w[4]+w[0]*(-exp(-(t-w[1])/w[3])+exp(-(t-w[1])/w[2]))/(-exp(-(w[1]+(w[2]*w[3]/(w[2]-w[3]))*ln(w[2]/w[3])-w[1])/w[3])+exp(-(w[1]+(w[2]*w[3]/(w[2]-w[3]))*ln(w[2]/w[3])-w[1])/w[2])))
	ENDIF
End
