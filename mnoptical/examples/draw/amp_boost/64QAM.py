import matplotlib.pyplot as plt
# import numpy as np


bpsk = [2.14558823825048e-08	,5.761659533770907e-10	,3.739417393705592e-13,	6.397451458013718e-14	]
qpsk = [5.356806832883978e-05	,8.38803712145313e-06	,1.9867147655802385e-07,	8.096794746262757e-08	]
eight_psk = [0.007474687756129199,	0.0032254733223997245,	0.0006017472374181713,	0.0004033276214277444]
# sixteen_psk = [0.21892490578024862	,0.2035277150301581,	0.17750642888253051	,0.17204917641629353]


gosnr = [15.006595,	18.524199,	25.703648,	27.441298]

plt.plot(gosnr, bpsk, marker='o',color='red', label="bpsk")
plt.plot(gosnr, qpsk, marker='o',color='blue', label="qpsk")
plt.plot(gosnr, eight_psk, marker='o',color='green', label="8psk")
# plt.plot(gosnr, sixteen_psk, marker='o',color='gray', label="16psk")

plt.legend(loc = 'upper right')

plt.xlabel('gosnr')
plt.title("gosnr and bit error rate in 64QAM")

plt.show()