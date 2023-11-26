
# FIXME: (AD) need to add all modulation formats with correct values

"""gOSNR threshold sensitivity defined in decibles (dB)"""
rx_thresholds = {
        '16QAM': 20,
        '64QAM': 10,
        '256QAM': 0
    }

"""bits per symbol"""
bps = {
    '16QAM': 4.0,
    '64QAM': 6.0,
    '256QAM':8.0
}

"""symbol rate"""
sr = {
    '16QAM': 32.0e9,
    '64QAM': 25.0e9,
    '256QAM':18.75e9
}