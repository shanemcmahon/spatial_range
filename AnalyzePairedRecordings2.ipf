

function CalcDelta(times, Ephys, BaseLinePeriod, ResponsePeriod, TimeShift)
	wave times,Ephys
	variable BaseLinePeriod,ResponsePeriod, TimeShift
	make /o /n=(numpnts(times)) DeltaV
	DeltaV = mean(Ephys,times[p]+TimeShift-BaseLinePeriod,times[p]+TimeShift) - mean(Ephys,times[p]+TimeShift,times[p]+TimeShift+ResponsePeriod)
End

macro DoCalcDelta(sTimes,sEphys,BaselinePeriod,ResponsePeriod,TimeShift)
	string sTimes, sEphys
	variable BaselinePeriod = 10e-3
	variable ResponsePeriod = 10e-3
  variable TimeShift = 0
	Prompt sTimes, "Event times name", popup, wavelist("*",";","")
	Prompt sEphys, "Ephys data name", popup, wavelist("*",";","")
	CalcDelta($sTimes, $sEphys,BaselinePeriod, ResponsePeriod, TimeShift)
endmacro

function ComparePDF(sTimes,sEphys,BaselinePeriod,ResponsePeriod,TimeShift)
	string sTimes, sEphys
	variable BaselinePeriod// = 2e-3
	variable ResponsePeriod// = 2e-3
  variable TimeShift// = 0
	Prompt sTimes, "Event times name", popup, wavelist("*",";","")
	Prompt sEphys, "Ephys data name", popup, wavelist("*",";","")
	CalcDelta($sTimes, $sEphys,BaselinePeriod, ResponsePeriod, 0)
	duplicate /o DeltaV DeltaV0
	CalcDelta($sTimes, $sEphys,BaselinePeriod, ResponsePeriod, TimeShift)
	// StatsKSTest /q  DeltaV0,DeltaV
	StatsTTest /q  DeltaV0,DeltaV
end

macro DoComparePDF(sTimes,sEphys,BaselinePeriod,ResponsePeriod,TimeShift)
	string sTimes, sEphys
	variable BaselinePeriod = 2e-3
	variable ResponsePeriod = 2e-3
  variable TimeShift = 0.05
	Prompt sTimes, "Event times name", popup, wavelist("*",";","")
	Prompt sEphys, "Ephys data name", popup, wavelist("*",";","")
	ComparePDF(sTimes,sEphys,BaselinePeriod,ResponsePeriod,TimeShift)
endmacro

Function BootNullKS(sTimes, sEphys,BaselinePeriod,ResponsePeriod,TimeShift,nBootReps)
string sTimes, sEphys
variable BaselinePeriod
variable ResponsePeriod
variable TimeShift
variable nBootReps
wave w_ksResults, w_statsttest
make /o /n=(nBootReps) wKSd
make /o /n=(nBootReps) wKSp
make /o /n=(nBootReps) wKSc
duplicate /o $sTimes, sTimesBoot
variable i
variable maxt = pnt2x($sEphys, numpnts($sEphys ))
variable mint = pnt2x($sEphys, 0)
for(i=0;i<nBootReps;i++)	// Initialize variables;continue test
sTimesBoot = mint + BaselinePeriod + ((maxt-mint)-(TimeShift + ResponsePeriod))*(0.5*(enoise(1)+1))
ComparePDF("sTimesBoot",sEphys,BaselinePeriod,ResponsePeriod,TimeShift)
// wKSd[i] = w_ksResults[4]
wKSd[i] = w_statsttest[8]
// wKSp[i] = w_ksResults[6]
wKSp[i] = w_statsttest[9]
endfor						// Execute body code until continue test is FALSE

ComparePDF(sTimes,sEphys,BaselinePeriod,ResponsePeriod,TimeShift)

end

macro DoBootNullKS(BaselinePeriod,ResponsePeriod,TimeShift,nBootReps,n,TimesCell,ResponseCell)

variable BaselinePeriod = 10e-3
variable ResponsePeriod = 10e-3
variable TimeShift = 0.05
variable nBootReps = 1000
variable n
string TimesCell = "cAt"
string ResponseCell = "cAr"
string sTimes, sEphys
sTimes = TimesCell+num2str(n)
sEphys = ResponseCell+num2str(n)
// Prompt sTimes, "Event times name", popup, wavelist("*",";","")
// Prompt sEphys, "Ephys data name", popup, wavelist("*",";","")
BootNullKS(sTimes, sEphys,BaselinePeriod,ResponsePeriod,TimeShift,nBootReps)
endmacro

dowindow /k d_statistic
dowindow /k p_val
dowindow /k dv_cdf

Display/K=0 /n=d_statistic wKSd
Display/K=0 /n=p_val wKSp
autopositionwindow /m=1 /r=d_statistic p_val

Display /n=dv_cdf /VERT DeltaV0,DeltaV
ModifyGraph /w=dv_cdf rgb(DeltaV0)=(0,0,0)
SetAxis /w=dv_cdf bottom -2,2
autopositionwindow /m=1 /r=p_val dv_cdf
Legend /w=dv_cdf /C/N=text0/F=0/M/H=5/A=MC
Legend /w=dv_cdf /C/N=text0/J/A=RB/X=0.00/Y=0.00

make /o/n=10 w_results = {w_ksResults[1],w_ksResults[4],w_ksResults[5],w_ksResults[6],wKSp(0.05),wKSp(0.5),wKSp(0.95),wKSd(0.05),wKSd(0.5),wKSd(0.95)}
make /o/t/n=10 w_results_names ={"n","d","c","p","p*5th","p*50th","p*95th","d*5th","d*50th","d*95th*"}

dowindow /k results
Edit/K=0 /n=results w_results_names,w_results
autopositionwindow /m=0 /r=d_statistic results
