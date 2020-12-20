from plot_imports import *
from math_imports import *
minor = 32; major = 36

def scaling_fun(x,a):
    z = a*(a+1)*x
    return a*(a+1)/gamma(a)/2.*pow(z,a/2.-1)*np.exp(-pow(z,0.5))

def mean_fun(t,a):
    return a*(a+1)/np.square(t) #s(0) = inf approximation
    #return 1./np.square(t/np.sqrt(a*(a+1)) + 1./np.sqrt(1)) # s(0) = 1

# first plot: scaling solution
ncols = 5
fig, axes, cbax = fig_with_cbar((32,6), ncols)
contents = np.loadtxt('frag_data_200.txt')
sel = np.array([0, 20, 60, 120, 200])
x = contents[0,:-2]; lam=contents[0,-2]; a=contents[0,-1]
#print(x)
c = contents[1:,1:-1][sel]
t = contents[1:,0][sel]
s = contents[1:,-1][sel]
print("delta function at x =",x[c[0]>0])
print("number of grid points:",len(x))
cmap = get_cmap(matplotlib.cm.get_cmap("plasma"),0,0.5)
cnorm = matplotlib.colors.Normalize(vmin=0, vmax=np.max(t))
for i,ax in enumerate(axes):
    ax.plot(x/s[i], s[i]*s[i]*c[i], color=cmap(cnorm(t[i])), alpha=0.5, lw=8, label='$\\phi_{num}$')
    ax.plot(x/s[-1], scaling_fun(x/s[-1],a), linestyle='dashed', color='k', alpha=0.8, lw=5, label='$\\phi_{ana}$')
    ax.set_xscale('log',basex=10)
    ax.set_yscale('log',basey=10)
    ax.set_ylim(7e-5,2e3)
    ax.set_xlim(2e-5,1e1)

    # set up individual figures too
    fig2 = plt.figure(figsize=(8,7))
    ax2 = fig2.add_subplot(111)
    ax2.plot(x/s[i], s[i]*s[i]*c[i], color=cmap(cnorm(40)), alpha=0.5, lw=8, label='$\\phi_{num}$')
    ax2.plot(x/s[-1], scaling_fun(x/s[-1],a), linestyle='dashed', color='k', alpha=0.8, lw=5, label='$\\phi_{ana}$')
    ax2.set_xscale('log',basex=10)
    ax2.set_yscale('log',basey=10)
    ax2.set_ylim(7e-8,2e4)
    ax2.set_xlim(8e-6,1e3)
    if i==0:
        set_legend(ax, 'lower left', minor)
        set_axis_labels(ax, "$\\xi$", "$\\phi(\\xi)$", major, xticks=[1e-3, 1e0, 1e3], yticks=[1e-6, 1e-3, 1e0, 1e3])
        set_legend(ax2, 'lower left', minor)
    else:
        set_axis_labels(ax, "$\\xi$", False, major, xticks=[1e-3, 1e0, 1e3])
    set_cbar(fig, cbax, cmap, cnorm, "$t$", [0, 5, 10, 15, 20], major)
    
    set_axis_labels(ax2, "$\\xi$", "$\\phi(\\xi)$", major, xticks=[1e-3, 1e0, 1e3], yticks=[1e-6, 1e-3, 1e0, 1e3])
    set_text(ax2, 2, 3e2, "$t="+str(int(t[i]))+"$", minor)
    fig2.subplots_adjust(hspace=0.35,wspace=0.05)
    fig2.tight_layout()
    fig2.savefig('vector_images/num_frag_'+str(i+1)+'.svg', format='svg', dpi=1200)

fig.tight_layout()
fig.savefig('vector_images/num_frag_tot.svg', format='svg', dpi=1200)

# second plot: s(t)
t = contents[1:,0]
s = contents[1:,-1]
fig = plt.figure(figsize=(8.5,8))
ax = fig.add_subplot(111)
tt = np.linspace(1e-3,500,1000)
ax.plot(t, s, color=cmap(cnorm(40)), alpha=0.5, lw=8, label='$s_{num}$')
ax.plot(tt, mean_fun(tt,a), linestyle='dashed', color='k', alpha=0.8, lw=5, label='$s_{ana}$')
ax.set_yscale('log')
ax.set_xscale('log')
#ax.set_xlim(-10,210)
ax.set_xlim(8e-1,500)
ax.set_ylim(8e-6,2)
set_axis_labels(ax, "$t$", "$s(t)$", major, xticks=[1e0, 1e1, 1e2], yticks=[1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0])
set_legend(ax, 'upper right', minor)
fig.tight_layout()
fig.savefig('vector_images/num_s_vs_t.svg', format='svg', dpi=1200)
plt.show()
