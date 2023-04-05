from collections import Counter

import numpy as np
import matplotlib.pyplot as plt

import xpart as xp
import xtrack as xt

fname_model = './generate_models/lhc_at_injection_obstacle.json'

collider = xt.Multiline.from_json(fname_model)
collider.build_trackers()

nemitt_x = 3.5e-6
nemitt_y = 3.5e-6
num_particles = 1000

bunch = xp.generate_matched_gaussian_bunch(line=collider.lhcb1,
    sigma_z=10e-2, nemitt_x=nemitt_x, nemitt_y=nemitt_y,
    num_particles=num_particles)

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

fig1.suptitle(
    r'$q_x$ = ' f'{tw["qx"]:.5f}' r' $q_y$ = ' f'{tw["qy"]:.5f}' '\n'
    r"$Q'_x$ = " f'{tw["dqx"]:.2f}' r" $Q'_y$ = " f'{tw["dqy"]:.2f}'
    r' $\gamma_{tr}$ = '  f'{1/np.sqrt(tw["momentum_compaction_factor"]):.2f}'
)

collider.lhcb1.track(bunch, num_turns=5)
print(f'Lost beam fraction: {np.sum(bunch.state <= 0)/num_particles*100:.2f}%')

# Get the lost particles
lost_particles = bunch.filter(bunch.state<=0)

# find the name of the element where we loose most particles
i_most_losses = Counter(lost_particles.at_element).most_common(1)[0][0]

print(f'The element where we loose most particles is called '
      f'`{collider.lhcb1.element_names[i_most_losses]}`.')

# Plot the lost particles and the obstacle
obstacle = collider.lhcb1['obstacle']
fig2 = plt.figure(2)
plt.plot(lost_particles.x, lost_particles.y, '.', markersize=1)
plt.plot(obstacle.x_vertices, obstacle.y_vertices, 'k')

s_obstacle = tw['s', 'obstacle']
element_close_to_obstacle = 'mb.c15r8.b1'

# Find vertical correctors close to the obstacle
tw.rows[s_obstacle-300.:s_obstacle+300.:'s', 'mcbv.*']

# Find out the knob that controls
collider.lhcb1.element_refs['mcbv.15r8.b1'].ksl[0]._expr

collider.lhcb1_co_ref.match(
    method='4d',
    ele_start='bpm.10r8.b1',
    ele_stop='bpm.18r8.b1',
    twiss_init=tw.get_twiss_init(at_element='bpm.10r8.b1'),
    vary=[
        xt.Vary(name='acbv11.r8b1', step=1e-10),
        xt.Vary(name='acbv13.r8b1', step=1e-10),
        xt.Vary(name='acbv15.r8b1', step=1e-10),
        xt.Vary(name='acbv17.r8b1', step=1e-10),
    ],
    targets=[
        # I want the vertical orbit to be at 10 mm the obstable with flat angle
        xt.Target('y', at=element_close_to_obstacle, value=10e-3, tol=1e-4, scale=1),
        xt.Target('py', at=element_close_to_obstacle, value=0, tol=1e-6, scale=1000),
        # I want the bump to be closed
        xt.Target('y', at='bpm.18r8.b1', value=tw['y', 'bpm.18r8.b1'],
                  tol=1e-6, scale=1),
        xt.Target('py', at='bpm.18r8.b1', value=tw['py', 'bpm.18r8.b1'],
                   tol=1e-7, scale=1000),
    ]
)

fig3 = plt.figure(3)
tw = collider.lhcb1.twiss()
spco = plt.subplot(2,1,1)
spco.plot(tw['s'], tw['y'], label='y')

plt.show()

