import numpy as np
import matplotlib.pyplot as plt

A1=10
EA1=17.45 
R1=287.05
T1=22
k1 = A1 ** (-EA1/(R1 * T1))
print(k1)
A2=101
EA2=17.45 
R2=287.05
T2=22
k2 = A2 ** (-EA2/(R2 * T2))
print(k2)
conc_A_0 = float(0.6)
time = np.arange(0.0, 60.0, 0.001)

conc_A = conc_A_0 * np.exp(-k1*time)
conc_C = conc_A_0 * (1 + (1/(k2-k1))*(k1 * np.exp(-k2 * time) - k2 * np.exp(-k1 * time)))

fig, ax = plt.subplots()

ax.plot(time, conc_A, label='concentration of Calcium Carbonate')
ax.plot(time, conc_C, label='concentration of Calcium Chloride')

ax.set_ylabel('concentration (mol/l)')
ax.set_xlabel('time (s)')
ax.legend()

plt.show()