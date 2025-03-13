import numpy as np
from utils.genStrlkFrame import gen_strlk_frame
from utils.getDiffEnc import get_diff_enc
from utils.ofdmSnrEstimator import ofdm_snr_estimator
from utils.getValidSC import get_valid_sc
from utils.fftAcqStrlk import fft_acq_strlk

# -------------------------------------- #
# --------- Starlink FRAME TEST -------- #
# -------------------------------------- #

def starlink_frame_test():
    # ----------- Test 1: No Noise
    N, Ng, Midx, Nsym, Nprovided = 1024, 32, 4, 300, 300
    s = {
        "SNRdB": np.nan,
        "Fsr": 240e6,
        "Fcr": 12075117187.5,
        "beta": 0,
        "type": "QAM",
        "data": np.random.randint(0, Midx, (N, Nprovided)),
    }
    
    y = gen_strlk_frame(s)
    y = np.concatenate((y, np.zeros(((N+Ng)*np.ceil(len(y)/(N+Ng)) - len(y),))))
    y = y.reshape((N+Ng, -1))
    y = y[Ng:, 2:302]
    Y = (1/np.sqrt(N)) * np.fft.fft(y)
    Y[:2, :] = np.nan
    Y[-2:, :] = np.nan
    
    if s["type"].upper() == "PSK":
        data_in = np.array([get_diff_enc(np.nan_to_num(np.exp(1j * np.pi * s["data"][:, i]))) for i in range(Nprovided)]).T
    elif s["type"].upper() == "QAM":
        data_in = np.array([get_diff_enc(np.nan_to_num(np.exp(1j * np.pi * s["data"][:, i] / np.max(s["data"][:, i])))) for i in range(Nprovided)]).T
    
    if Nsym - Nprovided > 0:
        Nreps, Nremain = divmod(Nsym, Nprovided)
        data_in = np.tile(data_in, (1, Nreps))
        data_in = np.concatenate((data_in, data_in[:, :Nremain]), axis=1)
    
    data_out = get_diff_enc(Y)
    match = np.sum(data_out == data_in, axis=0)
    str_result = "pass" if np.sum(match == N) == Nsym else "Failed"
    print(f'1. Starlink Frame Test no noise : {str_result:>30}')
    
    # Additional tests (Test 2-5) should follow the same pattern, adapting SNR handling, Doppler, and resampling.
    
if __name__ == "__main__":
    starlink_frame_test()