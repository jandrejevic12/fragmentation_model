from math_imports import *
from plot_imports import *
from setup_imports import *
from facet_utils import load_data, setup_dir

setup_dir(imdir)
len_inst = len_dist(a=0.0, name='len_dist')

def at_fun(s, t0):
    # start with an initial guess.
    tol = 1e-12
    a1 = 1.
    t = np.sqrt(a1*(a1+1)/s)
    a2 = np.sqrt(t/t0)
    while np.abs(a2-a1)>tol:
        a1 = a2
        t = np.sqrt(a1*(a1+1)/s)
        a2 = np.sqrt(t/t0)
    a = a2
    t = np.sqrt(a*(a+1)/s)
    return a, t

# Plot an example from manual data
Areas, Per, Lens = load_data('facet')
Shapes, Means, Times = load_data('props')
exp = 39
delta = 0.045
t = 2
fig, ax = plt.subplots(1,1,figsize=(9,8))
# data
lens = Lens[exp][t-1]/float(dim)
h = np.histogram(lens,bins2)[0]
hs = h/ds2/np.sum(h)
ax.scatter(sl2, hs, edgecolor=list(cmap(cnorm(delta)))[:-1]+[0.7], marker='o', s=200, lw=6, facecolor=list(cmap(cnorm(delta)))[:-1]+[0.4], zorder=10)
y = len_inst.pdf(sl2, Shapes[exp][t-1]+1, Means[exp][t-1], 0, 1)
ax.plot(sl2, y, color='royalblue', alpha=0.6, lw=6, label='best fit')
ax.set_xscale('log',basex=10)
ax.set_yscale('log',basey=10)
ax.set_ylim(7e-5,2e3)
ax.set_xlim(2e-5,1e1)
set_axis_labels(ax, "$r$", "$f_R(r)$", major, xticks=[1e-3, 1e-1, 1e1], yticks=[1e-3, 1e-1, 1e1, 1e3])
fig.tight_layout()
plt.savefig(imdir+'/len_plot_manual_'+str(exp)+'_'+str(t)+'.svg', format='svg', dpi=1200)
plt.show()
