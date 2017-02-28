macro AnalyzeGcamp()
  FindRespondingPixels()
  CalculateResponsePositions()
  CalculateDistanceToResponse()
  if (WaveExists($"root:AllPixelDistanceFrequency"))
 	   Make/N=10/O AllPixelDistance_Hist;DelayUpdate
		Histogram/B={0,0.5,10} AllPixelDistance,AllPixelDistance_Hist
		duplicate /o AllPixelDistance_Hist,root:temp
		root:AllPixelDistanceFrequency = root:AllPixelDistanceFrequency + root:temp
   else
   	Make/N=10/O AllPixelDistance_Hist;DelayUpdate
		Histogram/B={0,0.5,10} AllPixelDistance,AllPixelDistance_Hist
   	duplicate AllPixelDistance_Hist, root:AllPixelDistanceFrequency
   endif


EndMacro

function FindRespondingPixels()
string InputImageName
Prompt InputImageName, "Select Input Image", popup, WaveList("*", ";", "" )
DoPrompt /HELP="Select the image to use for counting responding pixels from the dropdown men"  "Select input image", InputImageName
wave InputImage = $InputImageName
variable i,j
make /o /n=(dimsize(InputImage,0),dimsize(InputImage,1)) RespondingPixels
RespondingPixels = 0
make /o /n=(DimSize(InputImage,2)) Tij
for(i=0;i<dimsize(InputImage,0);i++)
  for(j=0;j<DimSize(InputImage,1 );j++)	// Initialize variables;continue test
	 Tij = InputImage[i][j][p]
    RespondingPixels[i][j] = (mean(Tij,16,27) - mean(Tij,2,12)) > 4*(((2*Variance(Tij,2,12))^0.5)/(10^0.5))
  endfor						// Execute body code until continue test is FALSE
 endfor
end


function CalculateResponsePositions()

variable PixelSize
wave Parameters, RespondingPixels
Prompt PixelSize, "Enter pixel size (µm)"


if (WaveExists($"Parameters"))
 			// Execute if condition is TRUE
        // xPos = Parameters[0]
        // yPos = Parameters[1]
        PixelSize = Parameters[2]
      else
        PixelSize = 0.434216
        make /o /n=3 Parameters //<FALSE part>			// Execute if condition is FALSE
        Parameters[0] = 0.5
        Parameters[1] = 0.5
endif

DoPrompt /HELP="Enter uncaging position from Markpoints xml file. Allowed values are from 0-1" "Enter uncaging coordinates",PixelSize
Parameters[2] = PixelSize

SetScale /P x, 0,(PixelSize), RespondingPixels
SetScale /P y, 0,(PixelSize), RespondingPixels
duplicate /o RespondingPixels, RespondingPixelsX
duplicate /o RespondingPixels, RespondingPixelsY
RespondingPixelsX = x
RespondingPixelsY = y
Make /O /n=(0,2) RespondingPixelsCoordinates
Make /O /n=(0,2) AllPixelCoordinates
variable i,j
for(i=0;i<DimSize(RespondingPixels,0 );i++)	// Initialize variables;continue test
  for(j=0;j<DimSize(RespondingPixels,1 );j++)	// Initialize variables;continue test
     InsertPoints /m=0 0, 1, AllPixelCoordinates
     AllPixelCoordinates[0][0] = RespondingPixelsX[i][j]
     AllPixelCoordinates[0][1] = RespondingPixelsY[i][j]
    if (RespondingPixels[i][j] == 1)
     InsertPoints /m=0 0, 1, RespondingPixelsCoordinates
     RespondingPixelsCoordinates[0][0] = RespondingPixelsX[i][j]
     RespondingPixelsCoordinates[0][1] = RespondingPixelsY[i][j]
     else
        // no response, nothing to do
     endif

   endfor						// Execute body code until continue test is FALSE

endfor						// Execute body code until continue test is FALSE
end


function CalculateDistanceToResponse()
variable xPos, yPos, PixelSize
wave Parameters, RespondingPixels,RespondingPixelsX, RespondingPixelsY, RespondingPixelsCoordinates
wave AllPixelCoordinates
Prompt xPos,"Enter uncaging x position"
Prompt yPos,"Enter uncaging y position"
Prompt PixelSize, "Enter pixel size (µm)"

if (WaveExists($"Parameters"))
 			// Execute if condition is TRUE
        xPos = Parameters[0]
        yPos = Parameters[1]
        PixelSize = Parameters[2]
      else
        PixelSize = 0.434216
        xPos = 0.5
        yPos = 0.5
        make /o /n=3 Parameters //<FALSE part>			// Execute if condition is FALSE
endif

DoPrompt /HELP="Enter uncaging position from Markpoints xml file. Allowed values are from 0-1" "Enter uncaging coordinates",xPos,yPos,PixelSize
Parameters[0] = xPos
Parameters[1] = yPos
Parameters[2] = PixelSize
make /o /n=2 UncagingPosition
UncagingPosition[0] = PixelSize*xPos * Dimsize(RespondingPixels,0)
UncagingPosition[1] = PixelSize*yPos * Dimsize(RespondingPixels,1)
make /o /n=(DimSize(RespondingPixelsCoordinates,0 )) RespondingPixelsDistance
make /o /n=(DimSize(AllPixelCoordinates,0)) AllPixelDistance
RespondingPixelsDistance = sqrt((RespondingPixelsCoordinates[p][0]-UncagingPosition[0])^2 + (RespondingPixelsCoordinates[p][1]-UncagingPosition[1])^2)
AllPixelDistance = sqrt((AllPixelCoordinates[p][0]-UncagingPosition[0])^2 + (AllPixelCoordinates[p][1]-UncagingPosition[1])^2)

end
