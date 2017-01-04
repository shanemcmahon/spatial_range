edit//***********************************************************************
//***********************************************************************
//Original code from WS, modified by SM
//Implements Bootstrap test for difference in means on Studentized samples as described in Davison and Hinkley, p171 ex 4.19
//***********************************************************************
//***********************************************************************

//***********************************************************************

Function TValueCal(wave1, wave2) // fuction that calculate T value
    Wave wave1,wave2
    Variable Avg1,Var1,npnt1,Avg2,Var2,npnt2
    Variable Tvalue

    //gets stats from the 2 waves
    wavestats/Q wave1
    Avg1 = V_avg
    Var1 = ( V_sdev ) ^ 2
    npnt1 = V_npnts
    wavestats/Q wave2
    Avg2 = V_avg
    Var2 = ( V_sdev ) ^ 2
    npnt2 = V_npnts

    //calculate t
    Tvalue = (Avg1-Avg2) / sqrt ((Var1/npnt1)+(Var2/npnt2))
    Return Tvalue
END

//**************************************************************************

Macro BootstrapTtestFor2Samples (Oriwvname1, Oriwvname2, WvForT, resamplen)
    String WvForT,Oriwvname1, Oriwvname2
    Variable resamplen=1000
    Prompt Oriwvname1, "select the first sample wave", popup, WaveList("*",";","")//select wave1 from CDF
    Prompt Oriwvname2, "select the first sample wave", popup, WaveList("*",";","")//select wave2 from CDF
    Prompt WvForT, "name for the wave containing Tvalues"
    Prompt resamplen, "numbers of resampling trials:"

    silent 1
    PauseUpdate

    Variable BootstrapTrial=0
    Variable countTval=0
    Variable CalT=0
    Variable OriT=0
    Variable Pvalue=0
    Make/N=(resamplen)/O $WvForT



    Duplicate/O $Oriwvname1, ShiftWv1, BWave1
    wavestats/Q $Oriwvname1
    ShiftWv1=$Oriwvname1-V_avg         //create a new sample which has the mean at the zero while the SD is the same as the original sample
    BWave1=0 //wave for each bootstrap replication


    Duplicate/O $Oriwvname2, ShiftWv2, BWave2
    wavestats/Q $Oriwvname2
    ShiftWv2=$Oriwvname2-V_avg
    BWave2=0


    OriT = TValueCal ($Oriwvname1,$Oriwvname2) //calculate the t value for 2 original samples

       //************creating loop that draw resamples for 2 original samples seperately , and for each round, we calculate a t value for 2 re-samples
    do
        statsresample/N=(numpnts(ShiftWv1))/Q ShiftWv1 //resample command to generate a randomized resample of original wave with the same number of points
        BWave1=W_Resampled//create bootstrap replication for wave1
        statsresample/N=(numpnts(ShiftWv2))/Q ShiftWv2
        BWave2=W_Resampled//create bootstrap replication for wave2

        CalT= TValueCal(BWave1,BWave2)//calculate t value comparing BWave1 and BWave2 using TValueCal function


        $WvForT [BootstrapTrial] = CalT//write the t value for each bootstrap replication into container wave

        // If (abs(CalT)>abs(OriT)) // count how many times t  value is larger than the original t, 2-sided
        //     countTval+=1
        // Endif

        //If ((CalT)>(OriT)) // count how many times t  value is larger than the original t, 1-sided, µ1 > µ2
          If (abs(CalT)>abs(OriT)) // count how many times t  value is larger than the original t, 2-sided
            countTval+=1
        Endif

        BootstrapTrial+=1

    while (BootstrapTrial< resamplen)

    //*******************************

    Pvalue=countTval/resamplen // calculate the p value
    Killwaves ShiftWv1,ShiftWv2,BWave1,BWave2, W_Resampled
    Make/O/N=(resamplen/20)/O HistoBootstrapT
    Histogram/B=1 $WvForT,HistoBootstrapT
    Display HistoBootstrapT
    Cursor A, HistoBootstrapT, OriT
    Print "the waves have been tested are:", Oriwvname1,Oriwvname2
    Print "original t value is: ", OriT
    Print "bootstrap t test p value is:", Pvalue


Endmacro
