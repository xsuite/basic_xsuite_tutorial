import scipy.io as sio
import xtrack as xt

chamber_data = sio.loadmat('LHC_chm_ver.mat')
chamber = xt.LimitPolygon(
    x_vertices=chamber_data['Vx'][0],
    y_vertices=chamber_data['Vy'][0])

vx0 = chamber_data['Vx'][0].copy()
vy0 = chamber_data['Vy'][0].copy()

vx = vx0.copy()
vy = vy0.copy()

additional_points = [
    [-7e-3, -10e-3],
    [0e-3, -1.3e-3],
    [1e-2, 1e-3]
]

vx = list(vx[:-15])
vy = list(vy[:-15])

for pp in additional_points:
    vx.append(pp[0])
    vy.append(pp[1])

vx += list(vx0[-3:]) + [vx0[0]-1e-12]
vy += list(vy0[-3:]) + [vy0[0]-1e-12]

obstacle = xt.LimitPolygon(
    x_vertices=vx,
    y_vertices=vy)

fname_model = './lhc_at_injection.json'

collider = xt.Multiline.from_json(fname_model)

ulo_at = 'mb.c15r8.b1'
collider.lhcb1.insert_element(ulo_at, element=obstacle, name='obstacle')
collider.lhcb1.insert_element(ulo_at, element=chamber, name='chamber')

# Twiss, so that we can have the nice API to get the quads
collider.build_trackers()
tw = collider.lhcb1.twiss()
collider.lhcb1.discard_tracker()

# Put in the apertures around the quadrupoles near the obstacles
# Get quads for 300 meters around the obstacle each way
quads = tw.rows['obstacle%%-300':'obstacle%%300', 'mq.*'].cols['name'].rows[:]
# Insert apertures
for quad in quads:
    collider.lhcb1.insert_element(quad, element=chamber, name=f'chamber_for_{quad}')

collider.to_json('lhc_at_injection_obstacle.json')


