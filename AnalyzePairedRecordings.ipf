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
wave w_statsttest
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
  for(i=0;i<nEvents;i=i+1)	// Initialize variables;continue test
    // crosstalk = response - baseline
    CrossTalkRandom[i] = mean(wData2, (wLatencySim[i]+ResponseDelay-0.5*ResponsePeriod),(wLatencySim[i]+ResponseDelay+0.5*ResponsePeriod) ) - mean(wData2, (wLatencySim[i]-BaseLinePeriod),wLatencySim[i] )
  endfor						// Execute body code until continue test is FALSE
  StatsTTest /mean=0 /t=1 CrossTalkRandom
  TestStat1to2[j] = w_statsttest[5][0]
endfor						// Execute body code until continue test is FALSE




nEvents = numpnts(wLatency2)
make /o /n=(nEvents) CrossTalk2to1
for(i=0;i<nEvents;i=i+1)	// Initialize variables;continue test
// crosstalk = response - baseline
    CrossTalk2to1[i] = mean(wData1, (wLatency2[i]+ResponseDelay-0.5*ResponsePeriod),(wLatency2[i]+ResponseDelay+0.5*ResponsePeriod) ) - mean(wData2, (wLatency2[i]-BaseLinePeriod),wLatency2[i] )
endfor						// Execute body code until continue test is FALSE



end

menu "macros"
	"PairedAnaylsis"
end
