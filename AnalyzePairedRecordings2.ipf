

function CalcDelta(times, Ephys, BaseLinePeriod, ResponsePeriod, TimeShift)
	wave times,Ephys
	variable BaseLinePeriod,ResponsePeriod, TimeShift
	make /o /n=(numpnts(times)) DeltaV
	DeltaV = mean(Ephys,times[p]+TimeShift-BaseLinePeriod,times[p]+TimeShift) - mean(Ephys,times[p]+TimeShift,times[p]+TimeShift+ResponsePeriod)
End

macro DoCalcDelta(sTimes,sEphys,BaselinePeriod,ResponsePeriod,TimeShift)
	string sTimes, sEphys
	variable BaselinePeriod = 2e-3
	variable ResponsePeriod = 2e-3
  variable TimeShift = 0
	Prompt sTimes, "Event times name", popup, wavelist("*",";","")
	Prompt sEphys, "Ephys data name", popup, wavelist("*",";","")
	CalcDelta($sTimes, $sEphys,BaselinePeriod, ResponsePeriod, TimeShift)
endmacro
