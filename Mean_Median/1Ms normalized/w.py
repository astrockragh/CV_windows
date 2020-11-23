def sample_plot(sur, k, mean, sig, ax=None, legend=True, c='r',start=20, stop=-1, save=False, i=0, scale='log', top=None, bottom=1e-7):    
    start=20
    stop=-1
    c='r'
    h=cosmo.H0.value/100
    dz=[z_vals[k-1]-z_vals[k],z_vals[k]-z_vals[k+1]]
    norm=(cosmo.comoving_volume(z_vals[k]+dz[0]/2)-cosmo.comoving_volume(z_vals[k]-dz[1]/2))*area_dict[sur]
    norm=norm*h**3
    errcv=sig*mean
    errtot=np.sqrt(mean+(errcv)**2)
    poiserr=np.sqrt(mean)
    errcv, errtot, poiserr, mean=errcv/norm, errtot/norm, poiserr/norm, mean/norm
    print(errtot[k]/[errcv[k]])
    fig, ax=plt.subplots(figsize=(12,9))
    ax.plot(mass[start:stop], mean[k][start:stop], c, label=f'mass function for {sur}')
    ax.fill_between(mass[start:stop], mean[k][start:stop]+poiserr[k][start:stop],mean[k][start:stop]-poiserr[k][start:stop], color = c, alpha=0.4,
                label='Poisson uncertainty')
    ax.fill_between(mass[start:stop], mean[k][start:stop]+errcv[k][start:stop],mean[k][start:stop]-errcv[k][start:stop], color = c, alpha=0.2,
                label='Cosmic Variance')
    ax.fill_between(mass[start:stop], mean[k][start:stop]+ errtot[k][start:stop], mean[k][start:stop]- errtot[k][start:stop], color =c, alpha=0.1,
                label='Total statistical uncertainty')

    ax.yaxis.grid(color='gray', linestyle='dashed')
    # ax.hlines(1,8,11, 'k', linestyle='dashed')
    ax.set(xlabel=r'mass [$log(\frac{M_*}{M_\odot})$]', ylabel=r'n [$\frac{h^3}{Mpc^3}}$]',\
           title=f'mass function at z={np.round(z_vals[k],1)} for different survey designs', yscale=scale)
    ax.set_ylim(bottom, top)
    ax.legend(fontsize=10, loc='lower left')
    if save:
        plt.savefig(f'figs/masslumfunc/z_{np.round(z_vals[k],1)}_total.png')