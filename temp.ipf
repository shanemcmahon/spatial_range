macro do_save_results2()
save_results2()
endmacro

function save_results2()
wave rw2d_response, w_uncage_response,w2d_fake_pars
variable n_results
if(!waveexists(rw2d_response))
make /o /n=0 /t rw_uid

make /n=(numpnts(w_uncage_response),8) rw2d_response
make /o /n=8 rw_amp0_95
make /o /n=8 rw_amp0_90
make /o /n=8 rw_amp0_np_95
make /o /n=8 rw_amp0_np_90
endif
n_results = numpnts(rw_uid)
print n_results
setscale d,0,0,"pA",w_amplitude
make /o /n=(dimsize(w2d_fake_pars,0)) w_bs_amp0
w_bs_amp0[] = w2d_fake_pars[p][0]
sort w_bs_amp0, w_bs_amp0
setscale /i x,0,1,w_bs_amp0
setscale d,0,0,"pA",w_bs_amp0


make /o /n=(dimsize(w2d_fake_pars,0)) w_bs_amp0alt
w_bs_amp0alt[] = w2d_fake_pars[p][5]
sort w_bs_amp0alt, w_bs_amp0alt
setscale d,0,0,"pA",w_bs_amp0alt
setscale /i x,0,1,w_bs_amp0alt
print w_bs_amp0(0.1)
print w_bs_amp0(0.05)
print w_bs_amp0alt(0.1)
print w_bs_amp0alt(0.05)

// wave w2d_average_response
// variable /g n
// string wave_name
// Prompt wave_name,"uncaging response wave name",popup,wavelist("*",";","")
// doprompt "",wave_name
// wave w_response = $wave_name
// if(n == DimSize(w2d_average_response, 1))
// redimension /n=(-1,2*DimSize(w2d_average_response, 1)) w2d_average_response
// endif
//
// w2d_average_response[][n]=w_response[p]
// n=n+1
// movewave w2d_average_response root:results:
// movevariable n root:results:
end
