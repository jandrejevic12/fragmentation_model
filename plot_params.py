from math_imports import *
from plot_imports import *
from setup_imports import *
import pickle

cmap = matplotlib.cm.get_cmap("plasma")
cnorm = matplotlib.colors.Normalize(vmin=0, vmax=1)
markers = ['o','s','^','v']
Areas, Per, Lens = pickle.load(open('facet.p','rb'))
Shapes, Means, Times = pickle.load(open('props.p','rb'))
Areas_w, Per_w, Lens_w = pickle.load(open('facet_watershed_rev.p','rb'))
Shapes_w, Means_w, Times_w = pickle.load(open('props_watershed.p','rb'))

# ax1: l_meas vs (1-delta)*t/(a+1)
# ax2: l_model vs l_empir
# ax3: s vs t
fig1, axes1, cbax1 = fig_with_cbar((10,8),1)
fig2, axes2, cbax2 = fig_with_cbar((16,8),1)
fig3, axes3, cbax3 = fig_with_cbar((10,8),1)
ax1 = axes1[0]; ax2 = axes2[0]; ax3 = axes3[0]

x = np.linspace(0,60,100)
ax1.plot(x, A*x, linestyle='dashed', lw=7, color='k', alpha=0.5)
ax2.plot(2*x, 2*x, linestyle='dashed', lw=7, color='k', alpha=0.5)
t = np.logspace(0,2.5,100)
a = np.sqrt(t/tau)
ax3.plot(t, a*(a+1)/t**2, linestyle='dashed', lw=7, color='k', alpha=0.5)

for i,(exp,delta) in enumerate(zip(Exps,Deltas)):
    for t in range(len(Areas[exp])):
        if t<3 or exp in [45, 43, 41]: 
            per = (np.sum(Per[exp][t])-4*dim)/2./dim
            area = Areas[exp][t]/float(dim*dim)
            s = np.mean(area)
            tt = Times[exp][t]; a = Shapes[exp][t]
            ell = c1*(1-delta)*np.log(1+c2*ts[i][t]/delta)
            ax1.scatter(tt/(a+1)*(1-delta), per*(1-delta), edgecolor=list(cmap(cnorm(delta)))[:-1]+[0.9], marker=markers[t], s=200, lw=7, facecolor='w', zorder=10)
            ax2.scatter(ell, (1-delta)*A*tt/(a+1), edgecolor=list(cmap(cnorm(delta)))[:-1]+[0.9], marker=markers[t], s=200, lw=7, facecolor='w', zorder=10)
            ax3.scatter(tt, s, edgecolor=list(cmap(cnorm(delta)))[:-1]+[0.9], marker=markers[t], s=200, lw=7, facecolor='w', zorder=10)
       
ax3.set_xscale('log')
ax3.set_yscale('log')
set_axis_labels(ax1, "$(1-\\tilde{\\Delta})t / (a+1)$", "$\\ell_{meas.}$", major, xticks=[0, 20, 40, 60], yticks=[0, 40, 80, 120])
set_axis_labels(ax2, "$\\ell_{empir.}$", "$\\ell_{model}^{(t)}$", major, xticks=[0, 40, 80, 120], yticks=[0, 40, 80, 120])
set_axis_labels(ax3, "$t$", "$s(t)$", major, xticks=[1e0, 1e1, 1e2], yticks=[1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0])
set_cbar(fig1, cbax1, cmap, cnorm, "$\\tilde{\\Delta}$", [0.2,0.4,0.6,0.8], major)
set_cbar(fig2, cbax2, cmap, cnorm, "$\\tilde{\\Delta}$", [0.2,0.4,0.6,0.8], major)
set_cbar(fig3, cbax3, cmap, cnorm, "$\\tilde{\\Delta}$", [0.2,0.4,0.6,0.8], major)
fig1.tight_layout()
fig2.tight_layout()
fig3.tight_layout()
fig1.savefig('vector_images/l_meas_slope.svg', format='svg', dpi=1200)
fig2.savefig('vector_images/l_model_vs_l_empir.svg', format='svg', dpi=1200)
fig3.savefig('vector_images/s_vs_t.svg', format='svg', dpi=1200)
plt.show()
