import numpy as np
import healpy as hp
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import meander

def plot_zoom(scan, ra, dec, title, reso=3, var="pVal", range=[0, 6],cmap=None):
    """
    Make a zoomed-in skymap around a particular point (RA, decl) on the sky

    Parameters
    -----------
    scan: numpy array
        Healpix map of values to be plotted (usually a skymap)
    ra: float
        Right Ascension value at the center of the plot (best-fit RA or source RA)
    dec: float
        Declination value at the center of the plot (best-fit dec or source dec)
    title: str
        Plot title to use
    reso: float
        Resolution (arcmins), default 3
    cmap: matplotlib colormap or None
        colormap to use. if not set default is Seaborn "Blues"
    """
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
    """
    Adds a color bar to an existing healpy map

    Parameters
    -----------
    labels: float list
        list of points to be used (default [0., 2., 4., 6.])
    col_label: str
        label for colorbar (default IceCube Event Time)
    cmap: matplotlib colormap or None
        colormap to use. if not set default is "Blues"
    offset: int
        offset value for colorbar's label. default is -35
    """
    fig = plt.gcf()
    ax = fig.add_axes([0.95, 0.2, 0.03, 0.6])
    labels = labels
    cb = mpl.colorbar.ColorbarBase(ax, cmap="Blues" if cmap is None else cmap,
                        #norm=mpl.colors.Normalize(vmin=range[0], vmax=range[1]), 
                        orientation="vertical")

    cb.set_label(col_label, labelpad=offset, fontsize=18)
    cb.set_ticks([0., 1.])
    cb.set_ticklabels(labels)
    cb.update_ticks()

def plot_labels(src_dec, src_ra, reso):
    """
    Add labels to healpy zoom

    Parameters
    -----------
    src_dec: float
        Declination value at the center of the plot (best-fit dec or source dec)
    src_ra: float
        Right Ascension value at the center of the plot (best-fit RA or source RA)
    reso: float
        Resolution (arcmins)
    """
    fontsize = 20
    plt.text(-1*np.radians(1.75*reso),np.radians(0), r"%.1f$^{\circ}$"%(np.degrees(src_dec)),
             horizontalalignment='right',
             verticalalignment='center', fontsize=fontsize)
    plt.text(-1*np.radians(1.75*reso),np.radians(reso), r"%.1f$^{\circ}$"%(reso+np.degrees(src_dec)),
             horizontalalignment='right',
             verticalalignment='center', fontsize=fontsize)
    plt.text(-1*np.radians(1.75*reso),np.radians(-reso), r"%.1f$^{\circ}$"%(-reso+np.degrees(src_dec)),
             horizontalalignment='right',
             verticalalignment='center', fontsize=fontsize)
    plt.text(np.radians(0),np.radians(-1.75*reso), r"%.1f$^{\circ}$"%(np.degrees(src_ra)),
             horizontalalignment='center',
             verticalalignment='top', fontsize=fontsize)
    plt.text(np.radians(reso),np.radians(-1.75*reso), r"%.1f$^{\circ}$"%(-reso+np.degrees(src_ra)),
             horizontalalignment='center',
             verticalalignment='top', fontsize=fontsize)
    plt.text(np.radians(-reso),np.radians(-1.75*reso), r"%.1f$^{\circ}$"%(reso+np.degrees(src_ra)),
             horizontalalignment='center',
             verticalalignment='top', fontsize=fontsize)
    plt.text(-1*np.radians(2.35*reso), np.radians(0), r"declination", 
                ha='center', va='center', rotation=90, fontsize=fontsize)
    plt.text(np.radians(0), np.radians(-2.05*reso), r"right ascension", 
                ha='center', va='center', fontsize=fontsize)

def plot_events(dec, ra, sigmas, src_ra, src_dec, reso, sigma_scale=5., col = 'k', constant_sigma=False,
                    same_marker=False, energy_size=False, with_mark=True, with_dash=False,
                    label=''):
    """
    Adds events to a healpy zoom plot. Events are expected to be from self.llh.exp

    Parameters
    -----------
    dec: float array
        Array of declination values for each event
    ra: float array
        Array of Right Ascension values for each event
    sigmas: float array
        Angular error (circularized) to be plotted. Usually a 90% CL angular error
    src_ra: float
        Right Ascension value at the center of the plot (best-fit RA or source RA)
    src_dec: float
        Declination value at the center of the plot (best-fit dec or source dec)
    reso: float
        Resolution (arcmins)
    sigma_scale: float or None
        Value to rescale the sigma parameter (default 5.0).
        If None: no angular undertainty is plotted (sometimes used for very long time windows)
    col: str or str array
        color to use for each event plotted (default k=black)
    constant_sigma: bool
        Ignores sigma parameter and plots all markers with a size of 20.
    with_mark: bool
        Uses an x marker instead of o
    with_dash: bool
        Plot the angular error as a dashed contour.
        Usually used to indicated a removed event (e.g. alert event that triggered the analysis)
    same_marker, energy_size: bool
        Currently unused options.
    """
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
    """
    Load settings to be used as default plot settings.
    Includes Times New Roman font and size 12 font
    """
    # undo eventual matplotlibrc in user config
    mpl.rcdefaults()
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

def contour(ra, dec, sigma, nside):
    r""" Function for plotting contours on skymaps

    Parameters
    -----------
    ra: ndarray
        Array of ra for events 
    dec: ndarray
        Array of dec for events
    sigma: ndarray
        Array of sigma to make contours around events
    nside:
        nside of healpy map

    Returns
    --------
    Theta: array
        array of theta values of contour
    Phi: array
        array of phi values of contour 
    """
    dec = np.pi/2 - dec
    sigma = np.rad2deg(sigma)
    delta, step, bins = 0, 0, 0
    delta= sigma/180.0*np.pi
    step = 1./np.sin(delta)/20.
    bins = int(360./step)
    Theta = np.zeros(bins+1, dtype=np.double)
    Phi = np.zeros(bins+1, dtype=np.double)
    # define the contour
    for j in range(0,bins) :
        phi = j*step/180.*np.pi
        vx = np.cos(phi)*np.sin(ra)*np.sin(delta) + np.cos(ra)*(np.cos(delta)*np.sin(dec) + np.cos(dec)*np.sin(delta)*np.sin(phi))
        vy = np.cos(delta)*np.sin(dec)*np.sin(ra) + np.sin(delta)*(-np.cos(ra)*np.cos(phi) + np.cos(dec)*np.sin(ra)*np.sin(phi))
        vz = np.cos(dec)*np.cos(delta) - np.sin(dec)*np.sin(delta)*np.sin(phi)
        idx = hp.vec2pix(nside, vx, vy, vz)
        DEC, RA = hp.pix2ang(nside, idx)
        Theta[j] = DEC
        Phi[j] = RA
    Theta[bins] = Theta[0]
    Phi[bins] = Phi[0]

    return Theta, Phi

def plot_contours(proportions, samples):
    r""" Plot containment contour around desired level.
    E.g 90% containment of a PDF on a healpix map

    Parameters
    -----------
    proportions: list
        list of containment level to make contours for.
        E.g [0.68,0.9]
    samples: array
        array of values read in from healpix map
        E.g samples = hp.read_map(file)
        
    Returns
    --------
    theta_list: list
        List of arrays containing theta values for desired contours
    phi_list: list
        List of arrays containing phi values for desired contours
    """

    levels = []
    sorted_samples = list(reversed(list(sorted(samples))))
    nside = hp.pixelfunc.get_nside(samples)
    sample_points = np.array(hp.pix2ang(nside, np.arange(len(samples)))).T
    for proportion in proportions:
        level_index = (np.cumsum(sorted_samples) > proportion).tolist().index(True)
        level = (sorted_samples[level_index] +
                (sorted_samples[level_index+1] if level_index+1<len(samples) else 0))/2.0
        levels.append(level)
    contours_by_level = meander.spherical_contours(sample_points, samples, levels)

    theta_list = []; phi_list=[]
    for contours in contours_by_level:
        for contour in contours:
            theta, phi = contour.T
            phi[phi<0] += 2.0*np.pi
            theta_list.append(theta)
            phi_list.append(phi)

    return theta_list, phi_list

def make_public_zoom_skymap(skymap, events, ra, dec, with_contour=True, name='test'):
    """
    Make a zoomed skymap for public webpage (currently unused, under development)

    Parameters
    -----------
    skymap: array
        Healpy probability skymap (plotted as a colorscale)
    events: array
        Array of events to plot (expected from self.llh.exp)
    ra: float
        Right Ascension value at the center of the plot (best-fit RA or source RA)
    dec: float
        Declination value at the center of the plot (best-fit dec or source dec)
    with_contour: bool
        Plot skymap 90% contour (default True)
    name: str
        Event name, to save in filename
    """

    pdf_palette = sns.color_palette("Blues", 500)
    cmap = mpl.colors.ListedColormap(pdf_palette)

    plot_zoom(skymap, ra, dec, "", range = [0,10], reso=3., cmap = cmap)

    plot_events(events['dec'], events['ra'], events['sigma'], ra, dec, 2*6, sigma_scale=1.0,
                constant_sigma=False, same_marker=True, energy_size=True, col = 'black')
    
    if with_contour:
        probs = hp.pixelfunc.ud_grade(skymap, 64)
        probs = probs/np.sum(probs)
        ### plot 90% containment contour of PDF
        levels = [0.9]
        theta, phi = plot_contours(levels, probs)
        hp.projplot(theta[0], phi[0], linewidth=2., c='k')
        for i in range(1, len(theta)):
            hp.projplot(theta[i], phi[i], linewidth=2., c='k', label=None)

    #plt.scatter(0,0, marker='*', c = 'k', s = 130, label = label_str) 
    #plt.legend(loc = 2, ncol=1, mode = 'expand', fontsize = 18.5, framealpha = 0.95)
    plot_color_bar(range=[0,6], cmap=cmap, col_label='GW Map Probability', offset = -10,
                   labels=[f'{min(skymap):.1e}',f'{max(skymap):.1e}'])

    plt.savefig(f'./{name}_skymap_zoom_public.png', bbox_inches='tight', dpi=300)
    plt.close()