from cpymad.madx import Madx
import xtrack as xt

import xmask as xm
import xmask.lhc as xmlhc

# Import user-defined optics-specific tools
import optics_specific_tools_hlhc15 as ost

config_fname = 'config_collisions.yaml'
out_fname = 'lhc_at_collisions.json'

# Read config file
with open(config_fname, 'r') as fid:
    config = xm.yaml.load(fid)
config_mad_model = config['config_mad']

# Make mad environment
xm.make_mad_environment(links=config_mad_model['links'])

# Start mad
mad_b1b2 = Madx(command_log="mad_collider.log")
mad_b4 = Madx(command_log="mad_b4.log")

# Build sequences
ost.build_sequence(mad_b1b2, mylhcbeam=1)
ost.build_sequence(mad_b4, mylhcbeam=4)

# Apply optics (only for b1b2, b4 will be generated from b1b2)
ost.apply_optics(mad_b1b2, optics_file=config_mad_model['optics_file'])

# Build xsuite collider
collider = xmlhc.build_xsuite_collider(
    sequence_b1=mad_b1b2.sequence.lhcb1,
    sequence_b2=mad_b1b2.sequence.lhcb2,
    sequence_b4=mad_b4.sequence.lhcb2,
    beam_config=config_mad_model['beam_config'],
    enable_imperfections=config_mad_model['enable_imperfections'],
    enable_knob_synthesis=config_mad_model['enable_knob_synthesis'],
    pars_for_imperfections=config_mad_model['pars_for_imperfections'],
    ver_lhc_run=config_mad_model['ver_lhc_run'],
    ver_hllhc_optics=config_mad_model['ver_hllhc_optics'])

conf_knobs_and_tuning = config['config_knobs_and_tuning']
for kk, vv in conf_knobs_and_tuning['knob_settings'].items():
    collider.vars[kk] = vv

# Save to file
collider.to_json(out_fname)