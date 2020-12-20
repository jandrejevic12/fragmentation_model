from math_imports import *
from plot_imports import *
from setup_imports import *
import pickle

sc_inst = scaling_dist(a=0.0, name='scaling_dist')
Areas, Per, Lens = pickle.load(open('facet.p','rb'))
Shapes, Means, Times = pickle.load(open('props.p','rb'))
Areas_w, Per_w, Lens_w = pickle.load(open('facet_watershed_rev.p','rb'))
Shapes_w, Means_w, Times_w = pickle.load(open('props_watershed.p','rb'))

# fitted a vs t
fig, axes, cbax = fig_with_cbar((10,8), 1)
ax = axes[0]
markers = ['o','s','^','v']
tt = np.linspace(0,250)
ax.plot(tt, np.sqrt(tt/tau), linestyle='dashed', lw=7, color='k', alpha=0.5)

for i,(exp,delta) in enumerate(zip(Exps,Deltas)):
    print("Calculating exp",exp)
    for t in range(len(Areas[exp])): 
    	if t<3 or exp in [45, 43, 41]:   
	        area = Areas[exp][t]/float(dim*dim)
	        s = np.mean(area)
	        g, loc, scale = sc_inst.fit(area/s, 1, floc=0, fscale=1)
	        tt = np.sqrt(g*(g+1)/s)
	        ax.scatter(tt, g, edgecolor=list(cmap(cnorm(delta)))[:-1]+[0.9], marker=markers[t], s=200, lw=7, facecolor='w', zorder=10)
        
ax.set_xlim(-5,220)
ax.set_ylim(0,3.5)
set_axis_labels(ax, "$t$", "$a$", major, xticks=[0, 50, 100, 150, 200], yticks=[0, 1, 2, 3])
set_cbar(fig, cbax, cmap, cnorm, "$\\tilde{\\Delta}$", [0.2,0.4,0.6,0.8], major)
fig.tight_layout()
fig.savefig('vector_images/a_vs_t.svg', format='svg', dpi=1200)
plt.show()