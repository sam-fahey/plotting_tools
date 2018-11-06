#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

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
zen = np.cos((Dec + 90) * np.pi/180)
t_of_day = MJD - min(MJD)
azi_0 = t_of_day * 360. + t_of_day / 365.25 * 360
azi = (RA + azi_0) % 360

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
                    rownorm=False,
                    vlim=[None,None],
                    xlabel='azimuth',
                    ylabel='zenith angle'):
    
    fig, ax = plt.subplots()

    if type(bins) is int: n_zen_bins, n_azi_bins = bins, bins
    elif type(bins) is list and len(bins)==2: n_zen_bins, n_azi_bins = bins
    else: print("\'bins\' object must be type int or list"); exit

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
            print "Oh no!"
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

    plt.tight_layout()
    if show: plt.show()

plot_zen_v_azi(zen, azi, show=True, rownorm=True, bins=12)
