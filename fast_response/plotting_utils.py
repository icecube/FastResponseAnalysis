def plot_zoom(scan, ra, dec, title, reso=3, var="pVal", range=[0, 6],cmap=None):
    if cmap is None:
        pdf_palette = sns.color_palette("Blues", 500)
        cmap = mpl.colors.ListedColormap(pdf_palette)
    hp.gnomview(scan, rot=(np.degrees(ra), np.degrees(dec), 0),
                    cmap=cmap,
                    max=max(scan),
                    reso=reso,
                    title=title,
                    notext=True,
                    cbar=False
                    #unit=r""
                    )

    plt.plot(4.95/3.*reso*np.radians([-1, 1, 1, -1, -1]), 4.95/3.*reso*np.radians([1, 1, -1, -1, 1]), color="k", ls="-", lw=3)
    hp.graticule(verbose=False)
    plot_labels(dec, ra, reso)

def plot_color_bar(labels=[0.,2.,4.,6.], col_label=r"IceCube Event Time", range=[0,6], cmap=None, offset=-35):
    fig = plt.gcf()
    #ax = fig.add_axes([0.25, -0.03, 0.5, 0.03])
    ax = fig.add_axes([0.95, 0.2, 0.03, 0.6])
    labels = labels
    cb = mpl.colorbar.ColorbarBase(ax, cmap=ps_map if cmap is None else cmap,
                        #norm=mpl.colors.Normalize(vmin=range[0], vmax=range[1]), 
                        orientation="vertical")
    #cb.ax.minorticks_on()

    cb.set_label(col_label, labelpad=offset, fontsize=18)
    cb.set_ticks([0., 1.])
    cb.set_ticklabels(labels)
    cb.update_ticks()
    #cb.ax.get_xaxis().set_ticklabels(labels)

def plot_labels(src_dec, src_ra, reso):
    """Add labels to healpy zoom"""
    fontsize = 20
    plt.text(-1*np.radians(1.75*reso),np.radians(0), r"%.2f$^{\circ}$"%(np.degrees(src_dec)),
             horizontalalignment='right',
             verticalalignment='center', fontsize=fontsize)
    plt.text(-1*np.radians(1.75*reso),np.radians(reso), r"%.2f$^{\circ}$"%(reso+np.degrees(src_dec)),
             horizontalalignment='right',
             verticalalignment='center', fontsize=fontsize)
    plt.text(-1*np.radians(1.75*reso),np.radians(-reso), r"%.2f$^{\circ}$"%(-reso+np.degrees(src_dec)),
             horizontalalignment='right',
             verticalalignment='center', fontsize=fontsize)
    plt.text(np.radians(0),np.radians(-1.75*reso), r"%.2f$^{\circ}$"%(np.degrees(src_ra)),
             horizontalalignment='center',
             verticalalignment='top', fontsize=fontsize)
    plt.text(np.radians(reso),np.radians(-1.75*reso), r"%.2f$^{\circ}$"%(-reso+np.degrees(src_ra)),
             horizontalalignment='center',
             verticalalignment='top', fontsize=fontsize)
    plt.text(np.radians(-reso),np.radians(-1.75*reso), r"%.2f$^{\circ}$"%(reso+np.degrees(src_ra)),
             horizontalalignment='center',
             verticalalignment='top', fontsize=fontsize)
    plt.text(-1*np.radians(2.35*reso), np.radians(0), r"declination", 
                ha='center', va='center', rotation=90, fontsize=fontsize)
    plt.text(np.radians(0), np.radians(-2.05*reso), r"right ascension", 
                ha='center', va='center', fontsize=fontsize)

def plot_events(dec, ra, sigmas, src_ra, src_dec, reso, sigma_scale=5., col = 'k', constant_sigma=False,
                    same_marker=False, energy_size=False, with_mark=True, with_dash=False,
                    label=''):
    """Adds events to a healpy zoom, get events from llh."""
    cos_ev = np.cos(dec)
    tmp = np.cos(src_ra - ra) * np.cos(src_dec) * cos_ev + np.sin(src_dec) * np.sin(dec)
    dist = np.arccos(tmp)

    if sigma_scale is not None:
        sigma = np.degrees(sigmas)/sigma_scale
        sizes = 5200*sigma**2
        if constant_sigma:
            sizes = 20*np.ones_like(sizes)
        if with_dash:
            hp.projscatter(np.pi/2-dec, ra, marker='o', linewidth=2, 
                edgecolor=col, linestyle=':', facecolor="None", s=sizes, 
                alpha=1.0)
        else:
            hp.projscatter(np.pi/2-dec, ra, marker='o', linewidth=2, 
                edgecolor=col, facecolor="None", s=sizes, alpha=1.0)
    if with_mark:
        hp.projscatter(np.pi/2-dec, ra, marker='x', linewidth=2, 
            edgecolor=col, facecolor=col, s=60, alpha=1.0)

def load_plotting_settings():
    mpl.use('agg')
    mpl.rcParams['text.usetex'] = True
    try:
        mpl.rcParams['text.latex.unicode'] = True
    except:
        # new mpl doesn't like this rcParam
        pass
    mpl.rcParams['mathtext.rm'] = 'Times New Roman'
    mpl.rcParams['mathtext.it'] = 'Times New Roman:italic'
    mpl.rcParams['mathtext.bf'] = 'Times New Roman:bold'

    mpl.rc('font', family='serif', size=12)
    mpl.rcParams['xtick.labelsize'] = 16
    mpl.rcParams['ytick.labelsize'] = 16
    mpl.rcParams['xtick.major.size'] = 5
    mpl.rcParams['ytick.major.size'] = 5