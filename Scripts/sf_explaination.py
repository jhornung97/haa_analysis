import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern']
plt.rcParams['font.size'] = 12

def smooth_step(x, x0, k):
    return 1 / (1 + np.exp(-k*(x - x0)))

x = np.linspace(20, 30, 10000)
y_ideal = np.heaviside(x - 27, 1)
y_true = smooth_step(x, 27, 5)
y_sim = smooth_step(x, 26.8, 5)

plt.plot(x, y_ideal, color='navy', label='ideal')
plt.plot(x, y_true, color='dimgray', label='true')
plt.plot(x, y_sim, color='darkorange', label='simulated')

plt.xticks([27], [r"$p_T^{\mathrm{thres}}$"])

plt.xlabel(r'particle $p_T$ [GeV]')
plt.ylabel('trigger efficiency')

plt.gca().spines['left'].set_position(('outward',0))
plt.gca().spines['bottom'].set_position('zero')

plt.gca().spines['left'].set_linewidth(1.5)
plt.gca().spines['bottom'].set_linewidth(1.5)

plt.gca().tick_params(width=1.5)

plt.gca().spines['right'].set_color('none')
plt.gca().spines['top'].set_color('none')

plt.legend(frameon=False, loc='upper left')

plt.savefig('/web/jhornung/public_html/misc/trigger_efficiency.png', bbox_inches='tight')

plt.show()