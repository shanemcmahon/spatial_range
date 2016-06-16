#pragma rtGlobals=3		// Use modern global access method and strict wave access.

function SaveWaveList()
string sWaveList
prompt sWaveList,"Wave list string"
doprompt "",sWaveList
sWaveList = ReplaceString(",", sWaveList, ";")
print sWaveList
save /T/B sWaveList
end

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

function which(mask,[result])
wave mask
wave result
Extract /indx/o mask, w_which, mask
if(!paramisdefault(result))
duplicate /o w_which result
endif
end

function subset1d(w_in,[mask,IndexList,w_out])
wave w_in,mask,IndexList,w_out
if(paramisdefault(IndexList))
if(paramisdefault(mask))
	abort("Must specifiy mask or IndexList")
endif
make /o/free IndexList
which(mask,result=IndexList)
endif
make /o /n=(numpnts(IndexList)) w_subset
w_subset = w_in[IndexList[p]]
if(!paramisdefault(w_out))
	duplicate /o w_subset w_out
endif

end

function subset2d(w_in,dim,[mask,IndexList,w_out])
wave w_in,mask,IndexList,w_out
variable dim
if(paramisdefault(IndexList))
if(paramisdefault(mask))
	abort("Must specifiy mask or IndexList")
endif
if(numpnts(mask)!=dimsize(w_in,dim))
abort("Mask must have same length as dimsize(w_in,dim)")
endif
make /o/free IndexList
which(mask,result=IndexList)
endif
if(dim==0)
make /o /n=((numpnts(IndexList)),dimsize(w_in,1)) w_subset
w_subset = w_in[IndexList[p]][q]
if(!paramisdefault(w_out))
	duplicate /o w_subset w_out
endif
endif
if(dim==1)
make /o /n=(dimsize(w_in,1),(numpnts(IndexList))) w_subset
w_subset = w_in[p][IndexList[q]]
if(!paramisdefault(w_out))
	duplicate /o w_subset w_out
endif
endif
end

function subset(w_in,[dim,mask,IndexList,w_out])
wave w_in, mask, IndexList, w_out
variable dim
if(paramisdefault(IndexList))
if(paramisdefault(mask))
	abort("Must specifiy mask or IndexList")
endif
make /o/free IndexList
which(mask,result=IndexList)
endif
	switch(wavedims(w_in))
		case 1:
		if(paramisdefault(w_out))
		subset1d(w_in,IndexList=IndexList)
		return 0
		endif
		subset1d(w_in,IndexList=IndexList,w_out=w_out)
			return 0
		case 2:
		if(paramisdefault(dim))
		abort("must specifiy dim parameter for 2d wave")
		endif
		if(paramisdefault(w_out))
		subset2d(w_in,dim,IndexList=IndexList)
		return 0
		endif
		subset2d(w_in,dim,IndexList=IndexList,w_out=w_out)
		return 0

	endswitch
end

Function ColumnMins(InputMatrix)
	wave InputMatrix
	if(!DimSize(InputMatrix, 1 ))
		abort("ColumnMins Input must be 2d wave")
	endif
	variable NRows = DimSize(InputMatrix, 0 )
	variable NColumns = DimSize(InputMatrix, 1)
		make /o /n=(NColumns) wColumnMins
	variable i
		make /free /n=(NRows) wColmnI
		for(i = 0;i < NColumns;i += 1)	// for1
		wColmnI = InputMatrix[p][i]
		wColumnMins[i] = wavemin(wColmnI)
		endfor // for1

end//ColumnMins


function PlotCols(wIn)
wave wIn
variable ncols, i
display wIn[][0]
ncols = dimsize(wIn,1)
	for(i=0;i<ncols;i=i+1)	// Initialize variables;continue test
		appendtograph wIn[][i]
	endfor												// Execute body code until continue test is FALSE
end

menu "macros"
"SaveWaveList/5"
end