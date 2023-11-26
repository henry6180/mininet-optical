import matplotlib.pyplot as plt
# import numpy as np


# ber = {
#     'bpsk': {1.0429943153521592e-10, 1.043000556527073e-10},
#     'qpsk': {},
#     '8psk': {},
#     '16psk': {}
# }

bpsk = [1.043005477037504e-10,	5.02E-06,	0.000472093801,	0.06303676657]
qpsk = [3.501905731995346e-06,	0.00089573,	0.009689612576,	0.1396904528]
eight_psk = [0.0021758688943846497,	0.027279622	,0.08389874408	,0.3192363301]
sixteen_psk = [0.19691653239674933	,0.247034111,	0.2778365764,	0.3294061386]

gosnr = [20.192195	,9.751813,	5.466967,	1.170099]

plt.plot(gosnr, bpsk, marker='o',color='red', label="bpsk")
plt.plot(gosnr, qpsk, marker='o',color='blue', label="qpsk")
plt.plot(gosnr, eight_psk, marker='o',color='green', label="8psk")
plt.plot(gosnr, sixteen_psk, marker='o',color='gray', label="16psk")

plt.legend(loc = 'upper right')

plt.xlabel('gosnr')
plt.title("gosnr and bit error rate in 256QAM")

plt.show()