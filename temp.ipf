variable i_
i_ = 6
DeletePoints/M=1 i_,1, rw2d_fit_amplitude
DeletePoints/M=1 i_,1, rw2d_amplitude_0
DeletePoints/M=1 i_,1, rw2d_amplitude_0_np
DeletePoints/M=1 i_,1, rw2d_fit_amplitude_se
DeletePoints/M=1 i_,1, rw2d_fit_decay_time
DeletePoints/M=1 i_,1, rw2d_fit_onset_delay
DeletePoints/M=1 i_,1, rw2d_fit_rise_time
DeletePoints/M=1 i_,1, rw2d_fit_start_time
DeletePoints/M=1 i_,1, rw2d_fit_stop_time
DeletePoints/M=1 i_,1, rw2d_fit_t0
DeletePoints/M=1 i_,1, rw2d_fit_y0
DeletePoints/M=1 i_,1, rw2d_response
DeletePoints/M=1 i_,1, rw2d_uncage_time
DeletePoints/M=2 i_,1, rw3d_fits
DeletePoints/M=2 i_,1, rw3d_uncaging_response
DeletePoints i_,1, rwPockelsVoltage
DeletePoints i_,1, rw_uid



duplicate /o rw2d_fit_amplitude rw2dAmpNorm
edit rw2dampnorm
rowmeans(rw2dAmpNorm,naremove=1)
duplicate /o wRowMeans rwAmpNormRowMeans

rowmeans(rw2d_response)
