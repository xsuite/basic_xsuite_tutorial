import numpy as np
import matplotlib.pyplot as plt

import xpart as xp
import xtrack as xt

fname_model = './generate_models/lhc_at_collisions.json'

nemitt_x = 2.5e-6
nemitt_y = 2.5e-6

collider = xt.Multiline.from_json(fname_model)
collider.build_trackers()

# the two beamlines are accessible in collider.lhcb1, collider.lhcb2
collider.lhcb1.element_names
collider.lhcb1.elements

# Some examples
collider.lhcb1['ip1'] # marker element at one of the collision points (ATLAS experiment)

# A bending magnet
collider.lhcb1['mb.b18l3.b1..1']
collider.lhcb1['mb.b18l3.b1..1'].to_dict() # shows its properties

# A quadrupole
collider.lhcb1['mqwa.a4l3.b1..2']
collider.lhcb1['mqwa.a4l3.b1..2'].to_dict() # shows its properties

# Twiss (computes orbit, optics and other quantities of interest)
twb1 = collider.lhcb1.twiss()
twb2 = collider.lhcb2.twiss().reverse() # to have the two in the same reference frame

#### Measure tune

num_turns = 500

p_co = twb1.particle_on_co.copy()

collider.lhcb1.track(p_co, num_turns=num_turns, turn_by_turn_monitor=True)
mon_co = collider.lhcb1.record_last_track

p1 = twb1.particle_on_co.copy()
p1.x += 0.1e-3
p1.y += 0.2e-3
collider.lhcb1.track(p1, num_turns=num_turns, turn_by_turn_monitor=True)
mon1 = collider.lhcb1.record_last_track

tw = twb1
plt.close('all')

# Plot turn-by-turn data
fig100 = plt.figure(100, figsize=(6.4, 4.8*1.5))
spx = plt.subplot(2,1,1)
spy = plt.subplot(2,1,2, sharex=spx)
spx.plot(mon_co.x.T)
spx.plot(mon1.x.T)
spy.plot(mon_co.y.T)
spy.plot(mon1.y.T)
spy.set_xlabel('Turn')
spx.set_ylabel('x [m]')
spy.set_ylabel('y [m]')

# Plot the transverse spectrum
fig101 = plt.figure(101, figsize=(6.4, 4.8*1.5))
spx = plt.subplot(2,1,1)
spy = plt.subplot(2,1,2, sharex=spx)
freq_axis = np.fft.rfftfreq(n=num_turns)
spx.plot(freq_axis, np.abs(np.fft.rfft(mon1.x[0, :])))
spy.plot(freq_axis, np.abs(np.fft.rfft(mon1.y[0, :])))





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

# vertical lines at IPs with labels next to them at the top of the plot
for ip in ['ip1', 'ip2', 'ip5', 'ip8']:
    spbet.axvline(x=tw['s', ip], color='k', alpha=0.2, linestyle='--')
    spco.axvline(x=tw['s', ip], color='k', alpha=0.2, linestyle='--')


fig1.suptitle(
    r'$q_x$ = ' f'{tw["qx"]:.5f}' r' $q_y$ = ' f'{tw["qy"]:.5f}' '\n'
    r"$Q'_x$ = " f'{tw["dqx"]:.2f}' r" $Q'_y$ = " f'{tw["dqy"]:.2f}'
    r' $\gamma_{tr}$ = '  f'{1/np.sqrt(tw["momentum_compaction_factor"]):.2f}'
)


# Survey(s)
sv1 = collider.lhcb1.survey(element0='ip1')
sv2 = collider.lhcb2.survey(element0='ip1').reverse()

twb1_ir1 = twb1.rows[twb1['s', 'ip1']-200:twb1['s', 'ip1']+200:'s']
svb1_ir1 = sv1.rows[sv1['s', 'ip1']-200:sv1['s', 'ip1']+200:'s']

twb2_ir1 = twb2.rows[twb2['s', 'ip1']-200:twb2['s', 'ip1']+200:'s']
svb2_ir1 = sv2.rows[sv2['s', 'ip1']-200:sv2['s', 'ip1']+200:'s']

plt.figure(2)
plt.plot(svb1_ir1.Z, svb1_ir1.X, label='x', linestyle='--', color='b')
plt.plot(svb2_ir1.Z, svb2_ir1.X, label='y', linestyle='--', color='r')

plt.plot(svb1_ir1.Z, svb1_ir1.X + twb1_ir1.x, label='x', linestyle='-', color='b')
plt.plot(svb2_ir1.Z, svb2_ir1.X + twb2_ir1.x, label='y', linestyle='-', color='r')


plt.show()