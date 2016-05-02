#pragma rtGlobals=3		// Use modern global access method and strict wave access.

Function Kill_Wave_List(wave_list)
//
// Kill_Wave_List kills the waves listed in wave_list
//
String wave_list //A string containing a list of wave names
Variable n_items = ItemsInList(wave_list) //variable to hold the number of items in the wave list
Variable i // variable for iteration control
	for(i=0;i<n_items;i+=1)	// for1
		print StringFromList(i, wave_list)
		wave w = $StringFromList(i, wave_list)
		killwaves /z w
	endfor												//for1
end

Function RowMeans(InputMatrix,[NaRemove])
wave InputMatrix
variable NaRemove
if(!DimSize(InputMatrix, 1 ))
abort("RowMeans Input must be 2d wave")
endif
variable NColumns = DimSize(InputMatrix, 1 )
make /free/o /n=(NColumns) WTemp = 1/NColumns
MatrixOp/O WRowMeans= InputMatrix x WTemp

//if NaRemove is set then remove NaN values and recompute average
//print NaRemove
if(paramisdefault(NaRemove))
	NaRemove = 0
endif
if(NaRemove)
variable NRows = numpnts(WRowMeans)
variable i

	for(i = 0;i < NRows;i += 1)	// for1
	// first loop checks to see if the output wave value is NaN, if the entry is a valid number, then do nothing
	//else remove Na's and calculate new column mean
	//the nested for loop is much slower than matrix multiplication, so we only replace the column if it is NaN
		if(numtype(WRowMeans[i]))//if2
			make /free/o /n=(NColumns) WTemp
			WTemp = InputMatrix[i][p]
			variable j=0, NNaN=0, RowSum =0
			for(j=0;j<NColumns;j += 1)	//for2
				if (numtype(WTemp[j]))//if3
					NNaN += 1
				else//if3
					RowSum += WTemp[j]
				endif//if3
			endfor	//for2
				WRowMeans[i] = RowSum/(NColumns - NNaN)
		endif //if2
	endfor // for1
endif

end//RowMeans

Function ColumnMeans(InputMatrix,[NaRemove])
variable NaRemove
wave InputMatrix
if(!DimSize(InputMatrix, 1 ))
abort("ColumnMeans Input must be 2d wave")
endif
variable NRows = DimSize(InputMatrix, 0 )
make /free/o /n=(NRows) WTemp = 1/NRows
MatrixOp/O WColumnMeans=  WTemp^t x InputMatrix
//if NaRemove is set then remove NaN values and recompute average
//print NaRemove
if(paramisdefault(NaRemove))
	NaRemove = 0
endif
if(NaRemove)
variable NColumns = numpnts(WColumnMeans)
variable i

	for(i = 0;i < NColumns;i += 1)	// for1
	// first loop checks to see if the output wave value is NaN, if the entry is a valid number, then do nothing
	//else remove Na's and calculate new column mean
	//the nested for loop is much slower than matrix multiplication, so we only replace the column if it is NaN
		if(numtype(WColumnMeans[i]))//if2
			make /free/o /n=(NRows) WTemp
			WTemp = InputMatrix[p][i]
			variable j=0, NNaN=0, ColumnSum =0
			for(j=0;j<NRows;j += 1)	//for2
				if (numtype(WTemp[j]))//if3
					NNaN += 1
				else//if3
					ColumnSum += WTemp[j]
				endif//if3
			endfor	//for2
				WColumnMeans[i] = ColumnSum/(NRows - NNaN)
		endif //if2
	endfor // for1
endif
end//ColumnMeans
