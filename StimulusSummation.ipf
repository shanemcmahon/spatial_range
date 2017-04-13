macro StimSum()
duplicate /o :FitResults:ResponseWaves:w_Response0 PosA
duplicate /o :FitResults:ResponseWaves:w_Response1 PosB
PosB = PosB + :FitResults:ResponseWaves:w_Response2 + :FitResults:ResponseWaves:w_Response3
PosB = PosB/3
duplicate /o :FitResults:ResponseWaves:w_Response4 PosC
PosC = PosC + :FitResults:ResponseWaves:w_Response5 + :FitResults:ResponseWaves:w_Response6
PosC = PosC/3
duplicate /o :FitResults:ResponseWaves:w_Response7 PosABC_response
SetScale /p x,0, dimdelta(:UAdata:posNo_0,0),"s", PosA
SetScale /p x,0, dimdelta(:UAdata:posNo_0,0),"s", PosB
SetScale /p x,0, dimdelta(:UAdata:posNo_0,0),"s", PosC
SetScale /p x,0, dimdelta(:UAdata:posNo_0,0),"s", PosABC_response

Rotate x2pnt(PosA,2e-3), PosB,PosC
Rotate x2pnt(PosA,2e-3),PosC
if(!exists("vTemp"))
variable vTemp
endif
duplicate /o PosABC_response, PosABC_sum
PosABC_sum = PosA + PosB + PosC
vTemp = mean(PosABC_response,0,10e-3)
PosABC_response = PosABC_response - vTemp
vTemp =  mean(PosABC_sum,0,10e-3)
PosABC_sum = PosABC_sum - vTemp

dowindow /k Responses
display /n=Responses PosABC_response
AppendToGraph /c=(0,0,0) PosABC_sum
Legend/C/N=text0/F=0/M/H=5/A=MC
make /o  W_coef={-1-11,0.00098042,0.0092695,0.0039595,0,0.01}
FuncFit/H="000011" DiffTwoExp2 W_coef PosABC_response /D
duplicate /o w_coef ResponseCoef
FuncFit/H="000011" DiffTwoExp2 W_coef PosABC_sum /D
duplicate /o w_coef SumCoef
dowindow /k FitCoefs
edit /N=FitCoefs ResponseCoef,SumCoef
AutoPositionWindow /M=1 /r=Responses FitCoefs
endmacro
