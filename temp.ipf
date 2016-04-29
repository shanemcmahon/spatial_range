macro do_save_results2()
save_results2()
endmacro 

function save_results2()
wave w2d_average_response
variable /g n
string wave_name
Prompt wave_name,"uncaging response wave name",popup,wavelist("*",";","")
doprompt "",wave_name
wave w_response = $wave_name
if(n == DimSize(w2d_average_response, 1)) 
redimension /n=(-1,2*DimSize(w2d_average_response, 1)) w2d_average_response
endif

w2d_average_response[][n]=w_response[p]
n=n+1
movewave w2d_average_response root:results:
movevariable n root:results:
end


