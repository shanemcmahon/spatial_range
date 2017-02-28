Function GraphResponses()
wave ModelPredition = :FitResults:ModelPredictionWave2d
variable i

for(i=0;i<(dimsize(ModelPredition,0));i++)	// Initialize variables;continue test

PauseUpdate; Silent 1		// building window...
DoWindow /K $("Response" + num2str(i))
Display /N=$("Response" + num2str(i)) :FitResults:ResponseWave2d[i][*],:FitResults:ModelPredictionWave2d[i][*]
SetDrawEnv xcoord=bottom, ycoord=rel
DrawText 10E-3,0.1,"*"
ModifyGraph rgb(ModelPredictionWave2d)=(0,0,0)
ModifyGraph noLabel=2
ModifyGraph lblMargin=10
ModifyGraph axThick=0
ModifyGraph margin(left)=50,margin(bottom)=50,margin(top)=5,margin(right)=5,width=185
ModifyGraph height=80
 endfor						// Execute body code until continue test is FALSE

doWindow /K ResponseId
edit /n=ResponseId uid
ModifyTable width(uid)=200
end
