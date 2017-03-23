#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
	make /o /n=(750*1E5) c1
	make /o /n=(750*1E5) c2
	c1 = gnoise(1)
	c2 = gnoise(1)
	setscale /i x,0,750,c1
	setscale /i x,0,750,c2

function PairedAnaylsis()
	wave w_ksResults
	string sLatency1,sData1,sLatency2,sData2
	wave wLatency1,wData1,wLatency2,wData2
	variable nEvents
	Prompt sLatency1, "latency1 wave name" , popup, WaveList("*", ";", "" )
	Prompt sData1, "data1 wave name" , popup, WaveList("*", ";", "" )
	Prompt sLatency2, "latency2 wave name" , popup, WaveList("*", ";", "" )
	Prompt sData2, "data2 wave name" , popup, WaveList("*", ";", "" )
	DoPrompt  "", sLatency1, sData1, sLatency2, sData2
	wave wLatency1 = $sLatency1
	wave wLatency2 = $sLatency2
	wave wData1 = $sData1
	wave wData2 = $sData2

	variable BaseLinePeriod,ResponsePeriod,ResponseDelay,nBoot
	BaseLinePeriod = 1E-3
	ResponsePeriod = 1E-3
	ResponseDelay = 3E-3
	nBoot=10
	Prompt BaseLinePeriod, "Baseline Period"
	Prompt ResponsePeriod, "Response Period"
	Prompt ResponseDelay, " Response Delay"
	Prompt nBoot, "number of bootstrap estimates"
	DoPrompt "",BaseLinePeriod, ResponsePeriod, ResponseDelay, nBoot

	nEvents = numpnts(wLatency1)
	make /o /n=(nEvents) CrossTalk1to2


	variable i=0
	for(i=0;i<nEvents;i=i+1)	// Initialize variables;continue test
		// crosstalk = response - baseline
		CrossTalk1to2[i] = mean(wData2, (wLatency1[i]+ResponseDelay-0.5*ResponsePeriod),(wLatency1[i]+ResponseDelay+0.5*ResponsePeriod) ) - mean(wData2, (wLatency1[i]-BaseLinePeriod),wLatency1[i] )
	endfor						// Execute body code until continue test is FALSE

	duplicate /o  wLatency1, wLatencySim


	variable j = 0
	make /o /n=(nEvents) CrossTalkRandom
	make /o /n=(nBoot) TestStat1to2

	for(j=0;j<nBoot;j=j+1)	// Initialize variables;continue test
		print ((pnt2x(wData1,numpnts(wData1))-(ResponsePeriod + ResponseDelay)))
		print j
		wLatencySim = ((enoise(1)+1)/2) * ((pnt2x(wData1,numpnts(wData1))-(ResponsePeriod + ResponseDelay)))
    		//wLatencysim = wLatencySim - vResponseLength
		for(i=0;i<nEvents;i=i+1)	// Initialize variables;continue test
			// crosstalk = response - baseline
			CrossTalkRandom[i] = mean(wData2, (wLatencySim[i]+ResponseDelay-0.5*ResponsePeriod),(wLatencySim[i]+ResponseDelay+0.5*ResponsePeriod) ) - mean(wData2, (wLatencySim[i]-BaseLinePeriod),wLatencySim[i] )
		endfor						// Execute body code until continue test is FALSE
		duplicate /o CrossTalkRandom KStestData
		wavestats /q KStestData
		KStestData = KStestData - v_avg
		KStestData = KStestData/v_sdev
		StatsKSTest /q /CDFF=myNormalCDF KStestData
		TestStat1to2[j] = w_ksResults[2]
	endfor						// Execute body code until continue test is FALSE




	nEvents = numpnts(wLatency2)
	make /o /n=(nEvents) CrossTalk2to1
	for(i=0;i<nEvents;i=i+1)	// Initialize variables;continue test
		// crosstalk = response - baseline
		CrossTalk2to1[i] = mean(wData1, (wLatency2[i]+ResponseDelay-0.5*ResponsePeriod),(wLatency2[i]+ResponseDelay+0.5*ResponsePeriod) ) - mean(wData2, (wLatency2[i]-BaseLinePeriod),wLatency2[i] )
	endfor						// Execute body code until continue test is FALSE
end

function myNormalCDF(inx)
	variable inx
	//wave wMyCDFpar
	//return statsNormalCDF(inx,wMyCDFpar[0],wMyCDFpar[1])
	return statsNormalCDF(inx,0,1)
end

function MakeFakeResponse()
string sBaselineTemplate
Prompt sBaselineTemplate, "Select wave to use as template for baseline", popup, wavelist("*",";","")
DoPrompt  "",sBaselineTemplate
wave wBaseLineTemplate = $sBaselineTemplate

Variable vNevents, vRecordLength,vResponseLength
vNevents = 100; vRecordLength=10;vResponseLength=0.02
Prompt vNevents, "number of events"
Prompt vRecordLength, "Length of output (s)"
Prompt vResponseLength, "Length of synthetic response (s)"
DoPrompt "",vNevents,vRecordLength,vResponseLength

variable w0,w1,w2,w3,w4,w5
w0 = -10e-12
Prompt w0, "Response Amplitude"
w1 = 0
Prompt w1, "Response Delay"
w2=10e-3
Prompt w2, "Response Decay time"
w3=1e-3
Prompt w3, "Response Rise time"
w4=0
Prompt w4, "y0"
w5 = 0
Prompt w5, "Response Delay 2"
DoPrompt "",w0,w1,w2,w3,w4,w5
make /o /n=6 wCoef = {w0,w1,w2,w3,w4,w5}

variable vTemplateLength = dimdelta(wBaseLineTemplate,0)*dimsize(wBaseLineTemplate,0)

if (vRecordLength < vTemplateLength)
 Duplicate /o /r=(0,vRecordLength) wBaseLineTemplate, wResponseOut
 else
  make /o /n=(numpnts(wBaseLineTemplate),ceil(vRecordLength/vTemplateLength)) wResponseOut
  wResponseOut = wBaseLineTemplate[p]
  Redimension/N=(vRecordLength/dimdelta(wBaseLineTemplate,0)+1) wResponseOut
  setscale /i x,0,(vRecordLength),"s",wResponseOut
 endif

make /o /n=(vnEvents) wResponseTimes
make /o /n=(vResponseLength/dimdelta(wBaseLineTemplate,0)) wSingleResponse
setscale /i x,0,vResponseLength,"s",wSingleResponse
wSingleResponse = DiffTwoExp2(wCoef,x)
variable i, vtemp1, vtemp2, vtemp3
for(i=0;i<vnEvents;i=i)	// Initialize variables;continue test
    vtemp1 = (enoise(0.5) + 0.5)*vRecordLength - vResponseLength
    vtemp1 = max(vtemp1, vRecordLength-0.050)
    wResponseTimes[i] = vtemp1
    vtemp2 = x2pnt(wResponseOut,vtemp1)
    vtemp3 = numpnts(wSingleResponse )
    wResponseOut[vtemp2,(vtemp2+vtemp3)] = wResponseOut[p] + wSingleResponse[p-vtemp2-1]
 endfor						// Execute body code until continue test is FALSE


end



Function DiffTwoExp2(w,t) : FitFunc
	Wave w
	Variable t


	IF    ((t-(w[1]+w[5])) < 0 )
		return w[4]
	ELSE
		//Difference in exponential model from Schutter, Erik De. Computational modeling methods for neuroscientists. The MIT Press, 2009. Chapter 6
		variable Gsyn, Fnorm, TPeak,t0,td,tr,y0,UncageTime
		Gsyn = w[0]
		t0 = w[1]
		td = w[2]
		tr = w[3]
		y0 = w[4]
		UncageTime = w[5]
		TPeak = (UncageTime+t0) + (td*tr)/(td-tr)*ln(td/tr)
		fnorm = 1/(-exp(-(Tpeak-(UncageTime+t0))/tr) + exp(-(Tpeak-(UncageTime+t0))/td))
		return (y0 + Gsyn*fnorm*(exp(-(t-(UncageTime+t0))/td)-exp(-(t-(UncageTime+t0))/tr)))
		// return (w[4]+w[0]*(-exp(-(t-(w[1]+w[5]))/(w[3]))+exp(-(t-(w[1]+w[5]))/(w[2])))/(-exp(-((w[1])+((w[2])*(w[3])/((w[2])-(w[3])))*ln((w[2])/(w[3]))-(w[1]))/(w[3]))+exp(-((w[1])+((w[2])*(w[3])/((w[2])-(w[3])))*ln((w[2])/(w[3]))-(w[1]))/(w[2]))))
	ENDIF
End
menu "macros"
	"PairedAnaylsis"
  "MakeFakeResponse"
end
