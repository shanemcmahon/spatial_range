//make figures for gcamp mouse
#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
variable startX,stopX, startY
startX = 12
stopX = 256
startY = 100
variable ScanLinePeriod = 0.00117714556556446

make /n=(stopX-startX) Xcoord
Xcoord = x+startX
duplicate /o Xcoord, Ycoord
Ycoord = startY 

duplicate /o /r=[][2] FF, DF
DF = FF[p][2] - F0[p][2]
Display DF
newimage /n=ReferenceImage Reference
newimage /n=LineScanImage LineScan
make /o ProfileCenter, ProfileWidth = 0
edit ProfileCenter,ProfileWidth

duplicate XCoord, zeroes
zeroes = 0

AppendToGraph /T /W=ReferenceImage Ycoord vs Xcoord
ModifyGraph /W=ReferenceImage rgb(Ycoord)=(65535,65535,65535)
AppendToGraph /T /W=LineScanImage zeroes
ModifyGraph /W=LineScanImage rgb(zeroes)=(65535,65535,65535)


AppendToGraph/T /W=ReferenceImage Ycoord[(ProfileCenter[0]-ProfileWidth[0]/2),(ProfileCenter[0]+ProfileWidth[0]/2)] vs Xcoord[(ProfileCenter[0]-ProfileWidth[0]/2),(ProfileCenter[0]+ProfileWidth[0]/2)]
AppendToGraph /T /W=LineScanImage zeroes[(ProfileCenter[0]-ProfileWidth[0]/2),(ProfileCenter[0]+ProfileWidth[0]/2)]
ModifyGraph /W=LineScanImage lsize(zeroes#1)=4

Display /N = ProfilePlots MIP_Profile0[][2]
setscale /p x,0,ScanLinePeriod,("s"), MIP_Profile0

variable i_=2
AppendToGraph/T /W=ReferenceImage Ycoord[(ProfileCenter[i_]-ProfileWidth[i_]/2),(ProfileCenter[i_]+ProfileWidth[i_]/2)] vs Xcoord[(ProfileCenter[i_]-ProfileWidth[i_]/2),(ProfileCenter[i_]+ProfileWidth[i_]/2)]
AppendToGraph /T /W=LineScanImage zeroes[(ProfileCenter[i_]-ProfileWidth[i_]/2),(ProfileCenter[i_]+ProfileWidth[i_]/2)]
ModifyGraph /W=LineScanImage lsize( $("zeroes#"+num2str(i_+1)))=4
ChooseColor
ModifyGraph /W =LineScanImage rgb($("zeroes#"+num2str(i_+1)))=(v_red,v_green,v_blue)
ModifyGraph /W=ReferenceImage rgb($("Ycoord#"+num2str(i_+1)))=(v_red,v_green,v_blue)
AppendToGraph /W=ProfilePlots $("MIP_Profile"+num2str(i_))[][2]
ModifyGraph /W=ProfilePlots rgb($("MIP_Profile"+num2str(i_)))=(v_red,v_green,v_blue)
setscale /p x,0,ScanLinePeriod,("s"), $("MIP_Profile"+num2str(i_))