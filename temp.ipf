macro do_save_results2()
if(!exists("rw_uid"))
make /o /n=0 /t rw_uid
endif
InsertPoints numpnts(rw_uid), 1, rw_uid
rw_uid[numpnts(rw_uid)] = uid
save_results2()
endmacro

function save_results2()
wave rw2d_response, w_uncage_response,w2d_fake_pars,w_uncage_time
wave w_fit_start_time, w_fit_stop_time, w_amplitude, w_amplitude_se, w_t0
wave w_decay_time, w_rise_time, w_y0, w_onset_delay, w_amplitude_0
wave w_amplitude_0_alt, w2d_responses, w2d_fits
variable n_results

if(!waveexists(rw2d_response))
//make /o /n=0 /t rw_uid

make /n=(numpnts(w_uncage_response),8) rw2d_response
make /n=(numpnts(w_uncage_time),8) rw2d_uncage_time
make /n=(numpnts(w_fit_start_time),8) rw2d_fit_start_time
make /n=(numpnts(w_fit_stop_time),8) rw2d_fit_stop_time
make /n=(numpnts(w_amplitude),8) rw2d_fit_amplitude
make /n=(numpnts(w_amplitude_se),8) rw2d_fit_amplitude_se
make /n=(numpnts(w_t0),8) rw2d_fit_t0
make /n=(numpnts(w_decay_time),8) rw2d_fit_decay_time
make /n=(numpnts(w_rise_time),8) rw2d_fit_rise_time
make /n=(numpnts(w_y0),8) rw2d_fit_y0
make /n=(numpnts(w_onset_delay),8) rw2d_fit_onset_delay
make /n=(numpnts(w_amplitude_0),8) rw2d_amplitude_0
make /n=(numpnts(w_amplitude_0_alt),8) rw2d_amplitude_0_np
//make /n=(numpnts(w_amplitude_0_alt),8) rw2d_amplitude_0_alt
duplicate w2d_responses rw3d_uncaging_response
redimension /n=(-1,-1,8) rw3d_uncaging_response
duplicate w2d_fits rw3d_fits
redimension /n=(-1,-1,8) rw3d_fits
make /o /n=8 rw_amp0_95
make /o /n=8 rw_amp0_90
make /o /n=8 rw_amp0_np_95
make /o /n=8 rw_amp0_np_90
endif

n_results = numpnts(rw_uid)-1
//InsertPoints n_results, 1, rw_uid
//rw_uid[n_results] = uid

if(n_results >  dimsize( rw2d_response, 1)*3/4)
Redimension /N=(-1, 2*n_results) rw2d_response
Redimension /N=(-1, 2*n_results) rw2d_uncage_time
Redimension /N=(-1, 2*n_results) rw2d_fit_start_time
Redimension /N=(-1, 2*n_results) rw2d_fit_stop_time
Redimension /N=(-1, 2*n_results) rw2d_fit_amplitude
Redimension /N=(-1, 2*n_results) rw2d_fit_amplitude_se
Redimension /N=(-1, 2*n_results) rw2d_fit_t0
Redimension /N=(-1, 2*n_results) rw2d_fit_decay_time
Redimension /N=(-1, 2*n_results) rw2d_fit_rise_time
Redimension /N=(-1, 2*n_results) rw2d_fit_y0
Redimension /N=(-1, 2*n_results) rw2d_fit_onset_delay
Redimension /N=(-1, 2*n_results) rw2d_amplitude_0
Redimension /N=(-1, 2*n_results) rw2d_amplitude_0_np
Redimension /N=(-1,-1, 2*n_results) rw3d_uncaging_response
Redimension /N=(-1,-1, 2*n_results) rw3d_fits
Redimension /N=(2*n_results) rw_amp0_95
Redimension /N=(2*n_results) rw_amp0_90
Redimension /N=(2*n_results) rw_amp0_np_95
Redimension /N=(2*n_results) rw_amp0_np_90


endif


rw2d_response[][n_results] = w_uncage_response[p]
rw2d_uncage_time[][n_results] = w_uncage_time[p]
rw2d_fit_start_time[][n_results] = w_fit_start_time[p]
rw2d_fit_stop_time[][n_results] = w_fit_stop_time[p]
rw2d_fit_amplitude[][n_results] = w_amplitude[p]
rw2d_fit_amplitude_se[][n_results] = w_amplitude_se[p]
rw2d_fit_t0[][n_results] = w_t0[p]
rw2d_fit_decay_time[][n_results] = w_decay_time[p]
rw2d_fit_rise_time[][n_results] = w_rise_time[p]
rw2d_fit_y0[][n_results] = w_y0[p]
rw2d_fit_onset_delay[][n_results] = w_onset_delay[p]
rw2d_amplitude_0[][n_results] = w_amplitude_0[p]
rw2d_amplitude_0_np[][n_results] = w_amplitude_0_alt[p]
rw3d_uncaging_response[][][n_results]=w2d_responses[p][q]
rw3d_fits[][][n_results]=w2d_fits[p][q]


//setscale d,0,0,"pA",w_amplitude
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

print "initial amplitude 95th percentile", w_bs_amp0(0.05)
print "initial amplitude 90th percentile", w_bs_amp0(0.1)
print "nonparametric initial amplitude 95th percentile", w_bs_amp0alt(0.05)
print "nonparametric initial amplitude 90th percentile", w_bs_amp0alt(0.1)


return 1
rw_amp0_95[n_results] = w_bs_amp0(0.1)
rw_amp0_90[n_results] = w_bs_amp0(0.05)
rw_amp0_np_95[n_results] = w_bs_amp0alt(0.1)
rw_amp0_np_90[n_results] = w_bs_amp0alt(0.05)

end
