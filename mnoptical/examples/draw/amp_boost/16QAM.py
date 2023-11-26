import matplotlib.pyplot as plt
# import numpy as np

# h= [0.168, 0.170, 0.172, 0.174, 0.176, 0.178, 0.180, 0.182, 0.184]


sixteen_psk = [0.174114396	,0.171512	,0.168758,	0.175928,	0.184242]

# gosnr = [26.769557	,27.61282	,28.52941	,26.19989	,23.6825]
gosnr = [26.77	,27.61	,28.53	,26.2	,23.68]

# plt.plot(gosnr, bpsk, marker='o',color='red', label="bpsk")
# plt.plot(gosnr, qpsk, marker='o',color='blue', label="qpsk")
# plt.plot(gosnr, eight_psk, marker='o',color='green', label="8psk")
# plt.plot(gosnr, sixteen_psk, marker='o',color='gray', label="16psk")

# label = ["3dB", "5dB", "10dB", "15dB", "17dB"]
plt.bar(gosnr, sixteen_psk, tick_label=gosnr, width=0.3)

# plt.legend(loc = 'upper right')

plt.xlabel('gosnr(dB)')
plt.title("gosnr and bit error rate in 16QAM")

plt.show()