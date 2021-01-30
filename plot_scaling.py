from math_imports import *
from plot_imports import *
from setup_imports import *

setup_dir(imdir)
rw_inst = rw_gamma(a=-np.inf, name="rw_fold")

def w(delta):
    return R/np.sqrt(1-delta**2)

def dt_model(a,t,delta):
    n = t/(a+1)/2.
    return alpha*rw_inst.sf(w(delta),n,a+1,loc=0,scale=1)*(1-delta)/delta

def dl_model(a,t,delta):
    return (1-delta)*A/(a+1)*dt_model(a,t,delta)

def l_model(n,delta):
    c1 = 2/R**2
    c2 = alpha*R**2/np.sqrt(2*np.pi)
    tt = c1*(1-delta**2)*np.log(1+c2*n/delta/(1+delta))
    return (1-delta)*A/(np.sqrt(tt/tau)+1)*tt

def dl_empir(n,delta):
    return c1*c2*(1-delta)/delta/(1+c2*n/delta)

def l_empir(n,delta):
    return c1*(1-delta)*np.log(1+c2*n/delta)

cmap = matplotlib.cm.get_cmap("plasma")
cnorm = matplotlib.colors.Normalize(vmin=0, vmax=1)
markers = ['o','s','^','v']
Areas, Per, Lens = load_data('facet')
Shapes, Means, Times = load_data('props')
Areas_w, Per_w, Lens_w = load_data('facet_watershed')
Shapes_w, Means_w, Times_w = load_data('props_watershed')


fig1, axes1, cbax1 = fig_with_cbar((12,8),1)
fig2, axes2, cbax2 = fig_with_cbar((12,8),1)
fig3, axes3, cbax3 = fig_with_cbar((12,8),1)
fig4, axes4, cbax4 = fig_with_cbar((12,8),1)
ax1 = axes1[0]; ax2 = axes2[0]; ax3 = axes3[0]; ax4 = axes4[0]
x = np.linspace(0,120,100)
ax1.plot(x, x, linestyle='dashed', color='k', alpha=0.5, lw=5)
ax2.plot(x, x, linestyle='dashed', color='k', alpha=0.5, lw=5)
ax3.plot(x, x, linestyle='dashed', color='k', alpha=0.5, lw=5)
ax4.plot(x, x, linestyle='dashed', color='k', alpha=0.5, lw=5)

for i,(exp,delta) in enumerate(zip(Exps,Deltas)):
    # watershed
    for t in range(len(Areas_w[exp])):
        per = (np.sum(Per_w[exp][t])-4*dim)/2./dim
        tt = Times_w[exp][t]; a = Shapes_w[exp][t]
        if t<len(Areas_w[exp])-1:
            per_1 = (np.sum(Per_w[exp][t+1])-4*dim)/2./dim
            ax1.scatter((per_1-per)*(1-delta), dl_empir(t+1,delta), edgecolor=list(cmap(cnorm(delta)))[:-1]+[0.3], marker='o', s=150, lw=0, facecolor=list(cmap(cnorm(delta)))[:-1]+[0.3], zorder=10)
            ax3.scatter((per_1-per)*(1-delta), dl_model(a,tt,delta), edgecolor=list(cmap(cnorm(delta)))[:-1]+[0.3], marker='o', s=150, lw=0, facecolor=list(cmap(cnorm(delta)))[:-1]+[0.3], zorder=10)
        ax2.scatter(per*(1-delta), l_empir(t+1,delta), edgecolor=list(cmap(cnorm(delta)))[:-1]+[0.3], marker='o', s=150, lw=0, facecolor=list(cmap(cnorm(delta)))[:-1]+[0.3], zorder=10)
        ax4.scatter(per*(1-delta), l_model(t+1,delta), edgecolor=list(cmap(cnorm(delta)))[:-1]+[0.3], marker='o', s=150, lw=0, facecolor=list(cmap(cnorm(delta)))[:-1]+[0.3], zorder=10)
            
    # manual
    for t in range(len(Areas[exp])):
        if t<3 or exp in [45, 43, 41]: 
            per = (np.sum(Per[exp][t])-4*dim)/2./dim
            tt = Times[exp][t]; a = Shapes[exp][t]
            if t<len(Areas[exp])-1 and ts[i][t]+1 == ts[i][t+1]:
                per_1 = (np.sum(Per[exp][t+1])-4*dim)/2./dim
                #print(exp, ts[i][t], (per_1-per)*(1-delta), dl_empir(ts[i][t],delta), dl_model(a,tt,delta))
                ax1.scatter((per_1-per)*(1-delta), dl_empir(ts[i][t],delta), edgecolor=list(cmap(cnorm(delta)))[:-1]+[0.9], marker=markers[t], s=200, lw=6, facecolor='w', zorder=20)
                ax3.scatter((per_1-per)*(1-delta), dl_model(a,tt,delta), edgecolor=list(cmap(cnorm(delta)))[:-1]+[0.9], marker=markers[t], s=200, lw=6, facecolor='w', zorder=20)
            ax2.scatter(per*(1-delta), l_empir(ts[i][t],delta), edgecolor=list(cmap(cnorm(delta)))[:-1]+[0.9], marker=markers[t], s=200, lw=6, facecolor='w', zorder=20)
            ax4.scatter(per*(1-delta), l_model(ts[i][t],delta), edgecolor=list(cmap(cnorm(delta)))[:-1]+[0.9], marker=markers[t], s=200, lw=6, facecolor='w', zorder=20)
       
ax1.set_xscale('log'); ax1.set_yscale('log')
ax3.set_xscale('log'); ax3.set_yscale('log')
ax1.set_xlim(1e-1,1e2); ax1.set_ylim(1e-1,1e2)
#ax2.set_xlim(1e-1,150); ax2.set_ylim(1e-1,150)
ax3.set_xlim(1e-1,1e2); ax3.set_ylim(1e-1,1e2)
#ax4.set_xlim(1e-1,150); ax4.set_ylim(1e-1,150)
set_axis_labels(ax1, "$\\delta{\\ell}_{meas.}$", "$\\delta{\\ell}_{empir.}$",major)
set_axis_labels(ax3, "$\\delta{\\ell}_{meas.}$", "$\\delta{\\ell}_{model}$", major)
set_axis_labels(ax2, "$\\ell_{meas.}$", "$\\ell_{empir.}$", major)
set_axis_labels(ax4, "$\\ell_{meas.}$", "$\\ell_{model}^{(n)}$", major)
set_cbar(fig1, cbax1, cmap, cnorm, "$\\tilde{\\Delta}$", [0.2,0.4,0.6,0.8], major)
set_cbar(fig2, cbax2, cmap, cnorm, "$\\tilde{\\Delta}$", [0.2,0.4,0.6,0.8], major)
set_cbar(fig3, cbax3, cmap, cnorm, "$\\tilde{\\Delta}$", [0.2,0.4,0.6,0.8], major)
set_cbar(fig4, cbax4, cmap, cnorm, "$\\tilde{\\Delta}$", [0.2,0.4,0.6,0.8], major)
fig1.tight_layout()
fig2.tight_layout()
fig3.tight_layout()
fig4.tight_layout()
fig1.savefig(imdir+'/dl_empir_vs_dl_meas.svg', format='svg', dpi=1200)
fig2.savefig(imdir+'/l_empir_vs_l_meas.svg', format='svg', dpi=1200)
fig3.savefig(imdir+'/dl_model_vs_dl_meas.svg', format='svg', dpi=1200)
fig4.savefig(imdir+'/l_model_vs_l_meas.svg', format='svg', dpi=1200)
plt.show()
