#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

# Plotting code starts here
#   Required arrays:
#   - azi : azimuth values in degrees
#   - zen : cosine of zenith angle values
#
# Begin binning

def plot_zen_v_azi (zen, azi, bins=20, 
                    colormap='PuOr_r',
                    show=False,
                    savefile=None,
                    rownorm=True,
                    vlim=[None,None],
                    xlabel='azimuth',
                    ylabel=r'cos($\theta$)',
                    title=''):
    ''' Plot 2D zenith angle (y-axis) over azimuth (x-axis)
    :type    zen: 1d array
    :param   zen: zenith angle values (deg, rad, or cosine)

    :type    azi: 1d array
    :param   azi: azimuth values (deg or rad)

    :type   bins: int or list of two ints
    :param  bins: bins in each dimension; if list, number in azimuth by number in zenith
    '''
    fig, ax = plt.subplots()

    if type(bins) is int: n_zen_bins, n_azi_bins = bins, bins
    elif type(bins) is list and len(bins)==2: n_zen_bins, n_azi_bins = bins
    else: print("\'bins\' object must be type int or list"); exit

    if max(zen) > 4: zen = zen * np.pi / 180. # zen from deg to rad
    if min(zen) >= 0 and max(zen) > 1: zen = np.cos(zen) # zen from rad to cosine

    if max(azi) >= 10: # azi in degrees
        azi_max = 360
        azi_units = 'deg'
    elif max(azi) < 10: # azi in radians
        azi_max = 2*np.pi
        azi_units = 'rad'

    azi_bins = np.linspace(0, azi_max, n_azi_bins+1)
    zen_bins = np.linspace(-1, 1, n_zen_bins+1)
    H, x, y = np.histogram2d( azi, zen, bins=[azi_bins, zen_bins])
    X, Y = np.meshgrid(x,y)
    if rownorm:
        H_rownorm = H.copy()
        rowsums = np.array([sum([H[i][j] for i in range(len(H))]) for j in range(len(H[0]))])
        for col in range(len(H)): H_rownorm[col] = H[col]*1.*n_azi_bins / rowsums
        if vlim==[None,None]:
            vlim[0] = min([np.min(H_rownorm), 2-np.max(H_rownorm)])
            vlim[1] = 2-vlim[0]
        mesh = ax.pcolormesh(X, Y, H_rownorm.T, cmap=colormap, label='', vmin=vlim[0], vmax=vlim[1])
        cbar = fig.colorbar(mesh, pad=0., label='azimuth PDF w.r.t. mean')
    else:
        mesh = ax.pcolormesh(X, Y, H.T, cmap=colormap, label='', vmin=vlim[0], vmax=vlim[1])
        cbar = fig.colorbar(mesh, pad=0., label='counts')

    # Set labels
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if azi_units == 'rad':
        ax.set_xticks(np.linspace(0, 2*np.pi, 5))
        ax.set_xticklabels(['0', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'])
    elif azi_units == 'deg':
        ax.set_xticklabels(['%i$^{\circ}$'%a for a in ax.get_xticks()])
    plt.tight_layout()
    if show: plt.show()
    if savefile != None: plt.savefig(savefile)

def plot_polar_view (zen, azi, view,
                    bins=20, 
                    colormap='PuOr_r',
                    show=False,
                    savefile=None,
                    rownorm=True,
                    vlim=[None,None],
                    ylabel='azimuth PDF w.r.t. latitude mean',
                    title='',
                    rlabels=False):
    ''' Plot 2D zenith angle (y-axis) over azimuth (x-axis)
    :type     zen: 1d array
    :param    zen: zenith angle values (deg or rad)

    :type     azi: 1d array
    :param    azi: azimuth values (deg or rad)

    :type    view: string
    :param   view: 'south', 'north', or 'sky'; which polar view is plotted

    :type    bins: int or list of two ints
    :param   bins: bins in each dimension; if list, number in azimuth by number in zenith
    '''
    fig, ax = plt.subplots()

    if type(bins) is int: n_zen_bins, n_azi_bins = bins, bins
    elif type(bins) is list and len(bins)==2: n_zen_bins, n_azi_bins = bins
    else: print("\'bins\' object must be type int or list"); exit

    if max(zen) > 4: zen = zen * np.pi / 180. # zen from deg to rad
    if max(azi) > 7: azi = azi * np.pi / 180. # azi from deg to rad

    if view=='south':
        # convert zenith to latitude in southern hemisphere
        zen_mod = np.ma.masked_outside(zen, np.pi/2., 135.*np.pi/180)
        azi_mod = np.ma.compressed(np.ma.masked_array(azi, mask=zen_mod.mask))
        zen_mod = np.ma.compressed(zen_mod)
        lat = -2.*(135. - zen_mod * 180./np.pi)
        lat = np.sin((90 + lat) * np.pi/180.)
    if view=='north':
        # convert zenith to latitude in northern hemisphere
        zen_mod = np.ma.masked_outside(zen, 135.*np.pi/180, np.pi)
        azi_mod = np.ma.compressed(np.ma.masked_array(azi, mask=zen_mod.mask))
        zen_mod = np.ma.compressed(zen_mod)
        lat = 2.*(135. - zen_mod * 180./np.pi)
        lat = np.sin((90 + lat) * np.pi/180.)
    
    azi_bins = np.linspace(0, 2*np.pi, n_azi_bins+1)
    lat_bins = np.sin( np.linspace(0, np.pi/2, n_zen_bins+1))
    H, x, y = np.histogram2d( lat, azi_mod, bins=[lat_bins, azi_bins])
    if rownorm: H = [ H[i] / np.average(H[i]) for i in range(len(H)) ]
    ax = plt.subplot(111, polar=True)

    Theta, R = np.meshgrid(azi_bins, lat_bins)
    mesh = ax.pcolormesh(Theta, R, H, cmap=colormap, label='', vmin=vlim[0], vmax=vlim[1])
    cbar = fig.colorbar(mesh, pad=0.1)
    cbar.ax.set_ylabel(ylabel) 
    if not rlabels: ax.set_rticks([])

    if show: plt.show()
    if savefile != None: plt.savefig(savefile)

if __name__ == '__main__':
    # Load down-/up-going data
    MJD_d, E_d, err_d, RA_d, Dec_d = np.loadtxt('downgoing_events.txt', skiprows=1, unpack=True)
    MJD_u, E_u, err_u, RA_u, Dec_u = np.loadtxt('upgoing_events.txt', skiprows=1, unpack=True)
    # Combine down-/up-going for all-aky data arrays
    MJD = np.append(MJD_d, MJD_u)
    E = np.append(E_d, E_u)
    err = np.append(err_d, err_u)
    RA = np.append(RA_d, RA_u)
    Dec = np.append(Dec_d, Dec_u)
    
    # Work in zenith and azimuth
    zen = (Dec + 90) * np.pi/180
    t_of_day = MJD - min(MJD)
    azi_0 = t_of_day * 360. + t_of_day / 365.25 * 360
    azi = (RA + azi_0) % 360
    
    #plot_zen_v_azi(zen, azi, show=True, rownorm=True, bins=12)
    plot_polar_view(zen, azi, show=True, rownorm=True, bins=20, view='south')
    plot_polar_view(zen, azi, show=True, rownorm=True, bins=10, view='north')
