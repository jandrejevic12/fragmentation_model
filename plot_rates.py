from plot_imports import *
from math_imports import *
from setup_imports import *
from image_imports import *
from facet_utils import inflate, setup_dir

setup_dir(imdir)

def get_facets(exp, ts):
    Facets = {}
    Maps = {}
    Areas = {}

    for t in ts:
        mat_contents = sio.loadmat(data_prefix+"exp"+str(exp)+"/H"+str(t)+".mat")
        orig_data = mat_contents["H"+str(t)]
        border = np.isnan(orig_data).astype(np.int)
        if t == 1 and exp == 28:
            f = "labeled_data/Exp"+str(exp)+"/D"+str(t)+"_shifted.png"
            border[:-28,:] = border[28:,:]
        elif t ==1 and exp == 39:
            f = "labeled_data/Exp"+str(exp)+"/D"+str(t)+"_shifted.png"
            border[:-7,:] = border[7:,:]
        else:
            f = "labeled_data/Exp"+str(exp)+"/D"+str(t)+".png"
        im = imread(f, as_gray=True)
        labels = label(im.astype(np.int), connectivity=2)
        labels = inflate(labels, border)
        print(exp,t,'filter complete!')
        regions = regionprops(labels)
        inds = np.zeros((len(regions),2))
        areas = np.zeros(len(regions))
        cmap = -1*np.ones((N,N))
        facets = []
        for i,region in enumerate(regions):
            inds[i,:] = region.centroid
            areas[i] = region.area
            cmap[labels == region.label] = i
            facet = np.zeros((N,N))
            facet[labels == region.label] = 1
            facets += [facet]
        Facets[t] = facets
        Maps[t] = cmap
        Areas[t] = areas
    return Facets, Maps, Areas

def get_fragments(Facets, Maps, Areas, n, frac=0.5):
    Areas_rx = []
    Areas_fxy = []
    maxj = np.argsort(Areas[n])[::-1]
    for j in maxj:
        inds, counts = np.unique(Maps[n+1][(Maps[n+1]>=0)&(Facets[n][j]==1)], return_counts=True)
        inds = inds.astype(np.int)
        inds = inds[counts >= frac*Areas[n+1][inds]]
        # store area IF current facet fragmented
        if len(inds) > 1: # fragmentation occurred
            # exclude potential cases of all fragments larger than original
            if np.any((Areas[n+1][inds]/float(Areas[n][j])<1).astype(np.int)):
                Areas_rx += [Areas[n][j]/float(dim*dim)]
        # store ratio of child to parent facets (may have remained unfragmented)
        for area in Areas[n+1][inds]:
            if area/float(Areas[n][j])<=1: # exclude potential case of child larger than parent
                Areas_fxy += [area/float(Areas[n][j])]
    return Areas_rx, Areas_fxy

def log_pxy(x,a):
    return a*x + np.log(a+1)

for dt in [[1,2],[2,3]]:
    for i,exp in enumerate(Exps):
        delta = Deltas[i]
        Facets, Maps, Areas = get_facets(exp, dt)
        Areas_rx, Areas_fxy = get_fragments(Facets, Maps, Areas, dt[0])
        
        # rx plot
        fig, axes, cbax = fig_with_cbar((8.3,7), 1)
        ax = axes[0]
        h = np.histogram(Areas_rx,bins)[0] # counts of areas which fragmented per bin
        ht = np.histogram(Areas[dt[0]]/float(dim*dim),bins)[0] # previous counts in each bin - the number of "trials"
        hs = h/ht.astype(np.float) # fraction of facets which fragmented
        # variance in the proportion of successes for n trials is p*(1-p)/n
        var = hs*(1-hs)/ht
        ax.scatter(sl, hs, zorder=10, marker='o', lw=5, s=150, edgecolor=list(cmap(cnorm(delta)))[:-1]+[0.6], facecolor='w', label='$n='+str(ts[0])+'$')
        ax.errorbar(sl, hs, yerr=np.sqrt(var), zorder=9, color=list(cmap(cnorm(delta)))[:-1]+[0.6], lw=0, elinewidth=5, alpha=0.6)
        xx = np.logspace(-7,1,500)
        ax.set_xscale('log',basex=10)
        ax.set_yscale('log',basey=10)
        ax.set_ylim(3e-3,3e0)
        ax.set_xlim(3e-6,3e0)
        set_axis_labels(ax, "$x$", "$r(x)\\Delta{t}$", major, xticks=[1e-4, 1e-2, 1e0], yticks=[1e-2, 1e-1, 1e0])
        set_cbar(fig, cbax, cmap, cnorm, "$\\tilde{\\Delta}$", [0.2,0.4,0.6,0.8], major) 
        fig.tight_layout()
        plt.savefig(imdir+"/rx_"+str(exp)+"_"+str(dt[0])+".svg", format='svg', dpi=1200)

        # fxy plot
        fig, axes, cbax = fig_with_cbar((8.3,7), 1)
        ax = axes[0]
        h = np.histogram(Areas_fxy,bins)[0]
        hs = h/ds/np.sum(h)
        ax.scatter(sl, hs, zorder=10, marker='o', lw=3, s=150, edgecolor=list(cmap(cnorm(delta)))[:-1]+[0.6], facecolor=list(cmap(cnorm(delta)))[:-1]+[0.3], label='$n='+str(ts[1])+'$')

        # fit with curve fit (1 fitting param)
        coeff, cov = sp.optimize.curve_fit(log_pxy, np.log(sl[hs>0]), np.log(hs[hs>0]), p0=[-0.5])
        coeffs = [coeff[0], np.log(coeff[0]+1)]
        ax.plot(sl, sl**coeffs[0]*np.exp(coeffs[1]), linestyle='dashed', color='k', alpha=0.5, lw=5)

        ax.set_xscale('log',basex=10)
        ax.set_yscale('log',basey=10)
        ax.set_ylim(3e-2,3e4)
        ax.set_xlim(3e-6,7e0)
        set_text(ax, 7e-4, 4e3,"$\\beta=%.2f$"%(np.round(coeffs[0],2)), minor)
        set_axis_labels(ax, "$x/y$", "$\\rho(x/y)$", major, xticks=[1e-4, 1e-2, 1e0], yticks=[1e0, 1e2, 1e4])
        set_cbar(fig, cbax, cmap, cnorm, "$\\tilde{\\Delta}$", [0.2,0.4,0.6,0.8], major)
        fig.tight_layout()
        plt.savefig(imdir+"/fxy_"+str(exp)+"_"+str(dt[0])+".svg", format='svg', dpi=1200)
