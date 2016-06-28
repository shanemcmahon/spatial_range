duplicate /o wave1 ach_1
duplicate /o wave2 ach_3
duplicate /o /r=(0,w_fit_start_time[5]) ach_1 wResponseAvgOut
duplicate /o /r=(0,w_fit_start_time[5]) ach_3 wPowerOut

make /o /n= (dimsize(w2d_responses,1)) wTemp


wTemp = w2d_responses[1][p] + w2d_responses[2][p] + w2d_responses[4][p] + w2d_responses[7][p] + w2d_responses[11][p]
wTemp = wTemp/5
wResponseAvgOut[w_fit_start_pt[0], w_fit_start_pt[0]+numpnts(wTemp)-1] = wTemp[p-w_fit_start_pt[0]]

wTemp = w2d_responses[3][p] + w2d_responses[5][p] + w2d_responses[8][p] + w2d_responses[12][p]
wTemp = wTemp/4
wResponseAvgOut[w_fit_start_pt[1], w_fit_start_pt[1]+numpnts(wTemp)-1] = wTemp[p-w_fit_start_pt[1]]

wTemp = w2d_responses[6][p] + w2d_responses[9][p] + w2d_responses[13][p]
wTemp = wTemp/3
wResponseAvgOut[w_fit_start_pt[2], w_fit_start_pt[2]+numpnts(wTemp)-1] = wTemp[p-w_fit_start_pt[2]]

wTemp = w2d_responses[10][p] + w2d_responses[14][p]
wTemp = wTemp/2
wResponseAvgOut[w_fit_start_pt[3], w_fit_start_pt[3]+numpnts(wTemp)-1] = wTemp[p-w_fit_start_pt[3]]

wTemp = w2d_responses[0][p] + w2d_responses[15][p]
wTemp = wTemp/2
wResponseAvgOut[w_fit_start_pt[4], w_fit_start_pt[4]+numpnts(wTemp)-1] = wTemp[p-w_fit_start_pt[4]]

Save/T wResponseAvgOut,wPowerOut as "wResponseAvgOut++.itx"