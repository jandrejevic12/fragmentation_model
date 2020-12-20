from math_imports import *
from plot_imports import *
from setup_imports import *
import pickle

sc_inst = scaling_dist(a=0.0, name='scaling_dist')

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
Areas, Per, Lens = pickle.load(open('facet.p','rb'))
Shapes, Means, Times = pickle.load(open('props.p','rb'))
exp = 41
delta = 0.27
t = 4
fig, ax = plt.subplots(1,1,figsize=(9,8))
# data
area = Areas[exp][t-1]/float(dim*dim)
print("Mean check:",np.mean(area),Means[exp][t-1])
h = np.histogram(area/np.mean(area),bins)[0]
hs = h/ds/np.sum(h)
ax.scatter(sl, hs, edgecolor=list(cmap(cnorm(delta)))[:-1]+[0.7], marker='o', s=200, lw=6, facecolor=list(cmap(cnorm(delta)))[:-1]+[0.4], zorder=10)

# fit
g, loc, scale = sc_inst.fit(area/np.mean(area), 1, floc=0, fscale=1)
a, tt = at_fun(np.mean(area), tau)
n = 0.5*np.sqrt(Shapes[exp][t-1]/(Shapes[exp][t-1]+1)/np.mean(area))
th = 1./(2*n*(Shapes[exp][t-1]+1))
print("Shape check:",g, a, Shapes[exp][t-1])
print("t check:",tt, Times[exp][t-1], 1./th)
y = sc_inst.pdf(sl, Shapes[exp][t-1], 0, 1)
ax.plot(sl, y, color='royalblue', alpha=0.6, lw=7, label='best fit')
ax.set_xscale('log',basex=10)
ax.set_yscale('log',basey=10)
ax.set_ylim(7e-8,2e4)
ax.set_xlim(8e-6,1e3)
#set_legend(ax, 'lower left', minor)
set_axis_labels(ax, "$\\xi$", "$\\phi(\\xi)$", major, xticks=[1e-3, 1e0, 1e3], yticks=[1e-6, 1e-3, 1e0, 1e3])
fig.tight_layout()
plt.savefig('vector_images/dist_plot_manual_'+str(exp)+'_'+str(t)+'.svg', format='svg', dpi=1200)
plt.show()

# Plot an example from watershed data
Areas, Per, Lens = pickle.load(open('facet_watershed_rev.p','rb'))
Shapes, Means, Times = pickle.load(open('props_watershed.p','rb'))
exp = 41
delta = 0.27
t = 24
fig, ax = plt.subplots(1,1,figsize=(9,8))
# data
area = Areas[exp][t-1]/float(dim*dim)
print("Mean check:",np.mean(area),Means[exp][t-1])
h = np.histogram(area/np.mean(area),bins)[0]
hs = h/ds/np.sum(h)
ax.scatter(sl, hs, edgecolor=list(cmap(cnorm(delta)))[:-1]+[0.7], marker='o', s=200, lw=6, facecolor=list(cmap(cnorm(delta)))[:-1]+[0.4], zorder=10)

# fit
g, loc, scale = sc_inst.fit(area/np.mean(area), 1, floc=0, fscale=1)
a, tt = at_fun(np.mean(area), tau)
n = 0.5*np.sqrt(Shapes[exp][t-1]/(Shapes[exp][t-1]+1)/np.mean(area))
th = 1./(2*n*(Shapes[exp][t-1]+1))
print("Shape check:",g, a, Shapes[exp][t-1])
print("t check:",tt, Times[exp][t-1], 1./th)
y = sc_inst.pdf(sl, Shapes[exp][t-1], 0, 1)
y2 = sc_inst.pdf(sl, g, 0, 1)
ax.plot(sl, y, color='royalblue', alpha=0.6, lw=6, label='universal')
ax.plot(sl, y2, color='k', alpha=0.6, lw=5, linestyle='dashed', label='individual')
ax.set_xscale('log',basex=10)
ax.set_yscale('log',basey=10)
ax.set_ylim(7e-8,2e4)
ax.set_xlim(8e-6,1e3)
set_legend(ax, 'lower left', minor)
set_axis_labels(ax, "$\\xi$", "$\\phi(\\xi)$", major, xticks=[1e-3, 1e0, 1e3], yticks=[1e-6, 1e-3, 1e0, 1e3])
fig.tight_layout()
plt.savefig('vector_images/dist_plot_watershed_'+str(exp)+'_'+str(t)+'.svg', format='svg', dpi=1200)
plt.show()
