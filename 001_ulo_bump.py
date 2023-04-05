import numpy as np
import matplotlib.pyplot as plt

import xpart as xp
import xtrack as xt

fname_model = './generate_models/lhc_at_injection.json'

collider = xt.Multiline.from_json(fname_model)
collider.build_trackers()

nemitt_x = 3.5e-6
nemitt_y = 3.5e-6

bunch = xp.generate_matched_gaussian_bunch(line=collider.lhcb1,
    sigma_z=10e-2, nemitt_x=nemitt_x, nemitt_y=nemitt_y, num_particles=1000)

tw = collider.lhcb1.twiss()
sigmas = tw.get_betatron_sigmas(nemitt_x=nemitt_x, nemitt_y=nemitt_y)

plt.close('all')
fig1 = plt.figure(1, figsize=(6.4, 4.8*1.5))
spbet = plt.subplot(2,1,1)
spco = plt.subplot(2,1,2, sharex=spbet)
spbet.plot(tw['s'], tw['betx'])
spbet.plot(tw['s'], tw['bety'])
spbet.set_ylim(bottom=0)
spco.plot(tw['s'], tw['x'], label='x')
spco.plot(tw['s'], tw['y'], label='y')
spco.legend()
spbet.set_ylabel(r'$\beta_{x,y}$ [m]')
spco.set_ylabel(r'(Closed orbit)$_{x,y}$ [m]')

# vertical line at element mb.b15r8.b1
ulo_at = 'mb.c15r8.b1'
spbet.axvline(x=tw['s', ulo_at], color='k', alpha=0.2, linestyle='--')
spco.axvline(x=tw['s', ulo_at], color='k', alpha=0.2, linestyle='--')


fig1.suptitle(
    r'$q_x$ = ' f'{tw["qx"]:.5f}' r' $q_y$ = ' f'{tw["qy"]:.5f}' '\n'
    r"$Q'_x$ = " f'{tw["dqx"]:.2f}' r" $Q'_y$ = " f'{tw["dqy"]:.2f}'
    r' $\gamma_{tr}$ = '  f'{1/np.sqrt(tw["momentum_compaction_factor"]):.2f}'
)

import scipy.io as sio
chamber_data = sio.loadmat('LHC_chm_ver.mat')
obstacle = xt.LimitPolygon(
    x_vertices=chamber_data['Vx'][0],
    y_vertices=chamber_data['Vy'][0])

collider.lhcb1.track(bunch, ele_start=0, ele_stop=ulo_at)


fig_ulo = plt.figure(2)
ax_ulo = fig_ulo.add_subplot(111)
ax_ulo.plot(obstacle.x_vertices, obstacle.y_vertices)
ax_ulo.plot(bunch.x, bunch.y, '.', markersize=1)


plt.show()

