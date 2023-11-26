import matplotlib.pyplot as plt
# import numpy as np


sixteen_psk = [0.196477965,	0.246335,	0.276864,	0.327251]

gosnr = [20.30201	,9.86602,	5.567946	,1.268656]

# plt.plot(gosnr, sixteen_psk, marker='o',color='gray', label="16psk")

# plt.legend(loc = 'upper right')


plt.bar(gosnr, sixteen_psk, tick_label=gosnr, width=0.8)

plt.xlabel('gosnr(dB)')
plt.title("gosnr and bit error rate in 16QAM")

plt.show()