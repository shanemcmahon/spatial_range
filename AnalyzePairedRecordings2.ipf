

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
for(i=0;i<nBootReps;i++)	// Initialize variables;continue test
sTimesBoot = BaselinePeriod + (maxt-(TimeShift + ResponsePeriod))*(0.5*(enoise(1)+1))
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
