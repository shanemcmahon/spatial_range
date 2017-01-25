//transform waves to zero mean, normalized amplitude
posNo_0 = posNo_0 - y0w[0]
posNo_1 = posNo_1 - y0w[1]
posNo_2 = posNo_2 - y0w[2]
posNo_3 = posNo_3 - y0w[3]
posNo_4 = posNo_4 - y0w[4]
posNo_5 = posNo_5 - y0w[5]
posNo_6 = posNo_6 - y0w[6]
posNo_7 = posNo_7 - y0w[7]
posNo_8 = posNo_8 - y0w[8]
posNo_9 = posNo_9 - y0w[9]

posNo_0 = posNo_0/MaxAmp
posNo_1 = posNo_1/MaxAmp
posNo_2 = posNo_2/MaxAmp
posNo_3 = posNo_3/MaxAmp
posNo_4 = posNo_4/MaxAmp
posNo_5 = posNo_5/MaxAmp
posNo_6 = posNo_6/MaxAmp
posNo_7 = posNo_7/MaxAmp
posNo_8 = posNo_8/MaxAmp
posNo_9 = posNo_9/MaxAmp

//move waves to subfolder
MoveWave $StringFromList((itemsinlist(WaveList("*PosNo*", ";","" ))-1),WaveList("*PosNo*", ";","" )), $(":"+StringFromList(0,ListOfFolders,",")+":")
MoveWave $StringFromList((itemsinlist(WaveList("*PosNo*", ";","" ))-1),WaveList("*PosNo*", ";","" )), $(":"+StringFromList(1,ListOfFolders,",")+":")
MoveWave $StringFromList((itemsinlist(WaveList("*PosNo*", ";","" ))-1),WaveList("*PosNo*", ";","" )), $(":"+StringFromList(2,ListOfFolders,",")+":")
MoveWave $StringFromList((itemsinlist(WaveList("*PosNo*", ";","" ))-1),WaveList("*PosNo*", ";","" )), $(":"+StringFromList(3,ListOfFolders,",")+":")
MoveWave $StringFromList((itemsinlist(WaveList("*PosNo*", ";","" ))-1),WaveList("*PosNo*", ";","" )), $(":"+StringFromList(4,ListOfFolders,",")+":")
MoveWave $StringFromList((itemsinlist(WaveList("*PosNo*", ";","" ))-1),WaveList("*PosNo*", ";","" )), $(":"+StringFromList(5,ListOfFolders,",")+":")
MoveWave $StringFromList((itemsinlist(WaveList("*PosNo*", ";","" ))-1),WaveList("*PosNo*", ";","" )), $(":"+StringFromList(6,ListOfFolders,",")+":")
MoveWave $StringFromList((itemsinlist(WaveList("*PosNo*", ";","" ))-1),WaveList("*PosNo*", ";","" )), $(":"+StringFromList(7,ListOfFolders,",")+":")
MoveWave $StringFromList((itemsinlist(WaveList("*PosNo*", ";","" ))-1),WaveList("*PosNo*", ";","" )), $(":"+StringFromList(8,ListOfFolders,",")+":")
MoveWave $StringFromList((itemsinlist(WaveList("*PosNo*", ";","" ))-1),WaveList("*PosNo*", ";","" )), $(":"+StringFromList(9,ListOfFolders,",")+":")
