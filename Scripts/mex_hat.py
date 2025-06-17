import matplotlib.pyplot as plt
import numpy as np

def mexican_hat(x, mu_sq, lamb):
    return mu_sq/2 * x**2 + lamb/4 * x**4

plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 18

x = np.linspace(-1.75, 1.75, 100)
mex_hat = mexican_hat(x, -1, 1)
no_mex_hat = mexican_hat(x, 1, 1)

fig, ax = plt.subplots()

plt.plot(x, mex_hat, color='navy', label=r'$\mu^2 < 0$')
plt.plot(x, no_mex_hat, color='dimgray', label=r'$\mu^2 > 0$')
plt.xticks([-1, 1], ['$-v$', '$v$'])
plt.yticks([])

plt.ylim(np.min(mex_hat) - 0.1, np.max(mex_hat) + 0.1)

plt.gca().spines['left'].set_position('zero')
plt.gca().spines['bottom'].set_position('zero')

plt.gca().spines['left'].set_linewidth(1.5)
plt.gca().spines['bottom'].set_linewidth(1.5)

plt.gca().tick_params(width=1.5)

plt.gca().spines['right'].set_color('none')
plt.gca().spines['top'].set_color('none')

plt.xlabel('$\Phi$')
plt.ylabel('$V(\Phi)$', rotation=0)

plt.gca().xaxis.set_label_coords(1.05, 0.3)
plt.gca().yaxis.set_label_coords(0.5, 1.025)

plt.legend(frameon=False, bbox_to_anchor=(1.05, 1), loc='upper left')

plt.savefig('/web/jhornung/public_html/misc/both_potentials.png', bbox_inches='tight')

plt.show()