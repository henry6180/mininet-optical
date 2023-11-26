import matplotlib.pyplot as plt
# import numpy as np


bpsk = [2.1456008855078436e-08	,5.761858104325659e-10	,3.739417393705592e-13	,7.78195984523198e-14]
qpsk = [5.356823055516896e-05	,8.38818498970984e-06	,1.9867147655802385e-07,	8.944347146605729e-08]
eight_psk = [0.0074746980635896085,	0.0032254989898363043	,0.0006017472374181713,	0.0004216111208131014]
# sixteen_psk = [0.21892493264388402	,0.2035278523	,0.1775064289,	0.1726417007]

gosnr = [15.006589,	18.524166,	25.700854,	27.248791]

plt.plot(gosnr, bpsk, marker='o',color='red', label="bpsk")
plt.plot(gosnr, qpsk, marker='o',color='blue', label="qpsk")
plt.plot(gosnr, eight_psk, marker='o',color='green', label="8psk")
# plt.plot(gosnr, sixteen_psk, marker='o',color='gray', label="16psk")

plt.legend(loc = 'upper right')

plt.xlabel('gosnr')
plt.title("gosnr and bit error rate in 256QAM")

plt.show()