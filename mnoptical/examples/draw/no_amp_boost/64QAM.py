import matplotlib.pyplot as plt
# import numpy as np


# ber = {
#     'bpsk': {1.0429943153521592e-10, 1.043000556527073e-10},
#     'qpsk': {},
#     '8psk': {},
#     '16psk': {}
# }

bpsk = [1.043000556527073e-10,	5.02E-06,	0.0004720937073,	0.06303676101]
qpsk = [3.5018972952466627e-06,	0.00089573,	0.009689611556,	0.1396904457]
eight_psk = [0.002175866535770725,	0.027279618,	0.08389873984,	0.3192363215]
sixteen_psk = [0.1969165146655657,	0.247034108,	0.2778365748,	0.3294061373]

gosnr = [20.1922,9.751813,	5.466968,	1.170099]

plt.plot(gosnr, bpsk, marker='o',color='red', label="bpsk")
plt.plot(gosnr, qpsk, marker='o',color='blue', label="qpsk")
plt.plot(gosnr, eight_psk, marker='o',color='green', label="8psk")
plt.plot(gosnr, sixteen_psk, marker='o',color='gray', label="16psk")

plt.legend(loc = 'upper right')

plt.xlabel('gosnr')
plt.title("gosnr and bit error rate in 64QAM")

plt.show()