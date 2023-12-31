from mnoptical.dataplane import km, m, dB, dBm
from math import pi, sqrt
from scipy.special import erfc
import numpy as np
from mnoptical.units import db_to_abs, abs_to_db
from mnoptical.edfa_params import fibre_spectral_attenuation
from sys import argv

def calc(length=100, roadm_insertion_loss=17*dB,
         numAmp=2, boost_target_gain = 17*dB, ch=1):
    # Setting the parameter
    input_power=1e-3
    bw=32e09
    if ch==1:
        ch_freq = 191.35e12
    else:
        ch_freq = 191.40e12
    h = 6.62607015e-34
    Amp_gain=0
    span_loss = (list(fibre_spectral_attenuation['SMF']))[92-ch]
    if numAmp!=0:
        span_loss = span_loss*(length/numAmp*km)
        Amp_gain = length/numAmp*km*0.22
    else:
        span_loss = span_loss*(length*km)
    power = []
    # calculate the output power
    output_power = input_power*db_to_abs(-1*roadm_insertion_loss)
    output_power = output_power*db_to_abs(boost_target_gain)
    power+=[output_power]
    if numAmp!=0:
        for i in range(numAmp):
            output_power = output_power*db_to_abs(-1*span_loss)
            output_power = output_power*db_to_abs(Amp_gain)
            power+=[output_power]
    else:
        output_power = output_power*db_to_abs(-1*span_loss)
        power+=[output_power]
    
    # calculate the ase noise
    output_ase_noise = 0
    output_ase_noise =output_ase_noise*db_to_abs(boost_target_gain)+db_to_abs(5.5)*h*ch_freq*bw*db_to_abs(boost_target_gain)
    if numAmp!=0:
        for i in range(numAmp):
            output_ase_noise = output_ase_noise*db_to_abs(-1*span_loss)
            output_ase_noise = output_ase_noise*db_to_abs(Amp_gain)+db_to_abs(5.5)*h*ch_freq*bw*db_to_abs(Amp_gain)
    else:
        output_ase_noise = output_ase_noise*db_to_abs(-1*span_loss)

    # calculate the nli noise
    output_nli_noise = 0
    if numAmp!=0:
        for i in range(numAmp):
            output_nli_noise = output_nli_noise + gn_model(power[i],length/numAmp*km)
            output_nli_noise = output_nli_noise*db_to_abs(-1*span_loss)
            output_nli_noise = output_nli_noise*db_to_abs(Amp_gain)
    else:
        output_nli_noise = output_nli_noise + gn_model(power[0],length*km)
        output_nli_noise = output_nli_noise*db_to_abs(-1*span_loss)

    # calculate the last roadm part
    # since the roadm's carrier attenuation's calculation needs the input power, ase noise, nli noise
    carrier_attenuation = 0
    total_power = output_power+output_ase_noise+output_nli_noise
    carrier_attenuation = abs_to_db(total_power * 1e3) - (abs_to_db(input_power*1e03) - roadm_insertion_loss)
    if carrier_attenuation<0:
        carrier_attenuation=0.0
    output_power = output_power*db_to_abs(-1*carrier_attenuation)
    output_ase_noise = output_ase_noise*db_to_abs(-1*carrier_attenuation)
    output_nli_noise = output_nli_noise*db_to_abs(-1*carrier_attenuation)
    # calculate the last 0.001km span part
    last_span_loss = (list(fibre_spectral_attenuation['SMF']))[92-ch] *(0.001*km)
    output_nli_noise = output_nli_noise + gn_model(output_power,0.001*km)
    output_nli_noise = output_nli_noise*db_to_abs(-1*last_span_loss)
    output_power = output_power*db_to_abs(-1*last_span_loss)
    output_ase_noise = output_ase_noise*db_to_abs(-1*last_span_loss)
    
    # print("power=",output_power,", ","ase_noise=",output_ase_noise,", ","nli_noise=",output_nli_noise,". ")
    # print("OSNR=",abs_to_db(output_power/output_ase_noise),", ","gOSNR=",abs_to_db(output_power/(output_ase_noise+output_nli_noise)),". ")
    return [output_power, output_ase_noise, output_nli_noise, abs_to_db(output_power/output_ase_noise), abs_to_db(output_power/(output_ase_noise+output_nli_noise))]

def gn_model(power=1e-3, length=100):
    length = length * 1e03
    ref_wavelength=1550e-9
    dispersion = 1.67e-05
    non_linear_coefficient = 1.27 / 1e03
    attenuation_values = list(fibre_spectral_attenuation['SMF']) #ch1:0.21543, ch2:0.21532

    for i in range(0, len(attenuation_values)):
        attenuation_values[i] = attenuation_values[i] / 1e03
    fibre_attenuation = (attenuation_values)[::-1]
    alpha = fibre_attenuation / (20 * np.log10(np.e))
    beta2 = -(ref_wavelength ** 2) * abs(dispersion) / (2 * pi * 299792458.0 )
    gamma = non_linear_coefficient
    effective_length = (1 - np.exp(-2 * alpha * length)) / (2 * alpha)
    asymptotic_length = 1 / (2 * alpha)

    symbol_rate_cut = 32e09
    bw_cut = symbol_rate_cut
    pwr_cut = power
    g_cut = pwr_cut / bw_cut

    g_nli = 0
    symbol_rate_ch = 32e09
    bw_ch = symbol_rate_ch
    pwr_ch = power
    g_ch = pwr_ch / bw_ch
    psi = np.arcsinh(0.5 * np.pi ** 2 * asymptotic_length[0] * abs(beta2) * bw_cut ** 2)
    g_nli += g_ch ** 2 * g_cut * psi
    g_nli *= (16.0 / 27.0) * (gamma * effective_length[0]) ** 2 / (2 * np.pi * abs(beta2) * asymptotic_length[0])
    g_nli *= bw_cut
    return g_nli

def get_ber(gosnr):
    if gosnr>0:
        return (3/8)*erfc(sqrt(gosnr/10))
    return 1

if __name__ == '__main__':

    length=100
    roadm_insertion_loss=17*dB
    ber=1
    if len(argv)>=2: length = int(argv[1])
    if len(argv)>=3: roadm_insertion_loss = int(argv[2])*dB
    if len(argv)>=4: ber = float(argv[3])
    maxD=0.1 #max = maxD amps/km
    minB = int(argv[2])-5
    maxB = int(argv[2])+5
    maxA = int(length*maxD)
    if maxA<=10: maxA=10

    fo = open("result2.txt","w")
    for boost_target_gain in range(minB, maxB+1):
        for numAmp in range(1,maxA+1):
            t2_gosnr = (calc(length=length, roadm_insertion_loss=roadm_insertion_loss,
                             numAmp=numAmp, boost_target_gain=boost_target_gain*dB, ch=1))[4]
            t2_ber = get_ber(t2_gosnr)
            if t2_ber>ber:
                continue
            t1_gosnr = (calc(length=length, roadm_insertion_loss=roadm_insertion_loss,
                             numAmp=numAmp, boost_target_gain=boost_target_gain*dB, ch=2))[4]
            t1_ber = get_ber(t1_gosnr)
            if t2_ber>ber:
                continue
            fo.write("%d %d %.4f %.4f %.4f %.4f\n" % (boost_target_gain,numAmp,t1_gosnr,t2_gosnr,t1_ber,t2_ber))
            # fo.write(f'{boost_target_gain:<17d} {numAmp:<6d} ')
            # fo.write(f'{(t1_gosnr):8.4f} {(t2_gosnr):8.4f} {(t1_ber):6.4f} {(t2_ber):6.4f}\n')
    fo.close()