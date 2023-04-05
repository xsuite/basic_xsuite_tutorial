from xmask.lhc import install_errors_placeholders_hllhc

def build_sequence(mad, mylhcbeam, **kwargs):

    # Select beam
    mad.input(f'mylhcbeam = {mylhcbeam}')

    mad.input(

      f'''
      ! Get the toolkit
      call,file=
        "acc-models-lhc/toolkit/macro.madx";
      '''
      '''
      ! Build sequence
      option, -echo,-warn,-info;
      if (mylhcbeam==4){
        call,file="acc-models-lhc/lhc/lhcb4.seq";
      } else {
        call,file="acc-models-lhc/lhc/lhc.seq";
      };
      option, -echo, warn,-info;
      '''
      f'''
      !Install HL-LHC
      call, file=
        "acc-models-lhc/hllhc_sequence.madx";
      '''
      '''
      ! Slice nominal sequence
      exec, myslice;
      ''')

    install_errors_placeholders_hllhc(mad)

    mad.input(
      '''
      !Cycling w.r.t. to IP3 (mandatory to find closed orbit in collision in the presence of errors)
      if (mylhcbeam<3){
        seqedit, sequence=lhcb1; flatten; cycle, start=IP3; flatten; endedit;
      };
      seqedit, sequence=lhcb2; flatten; cycle, start=IP3; flatten; endedit;

      ! Install crab cavities (they are off)
      call, file='acc-models-lhc/toolkit/enable_crabcavities.madx';
      on_crab1 = 0;
      on_crab5 = 0;

      ! Set twiss formats for MAD-X parts (macro from opt. toolkit)
      exec, twiss_opt;


      '''
    )


def apply_optics(mad, optics_file):
    mad.call(optics_file)
    # A knob redefinition
    mad.input('on_alice := on_alice_normalized * 7000./nrj;')
    mad.input('on_lhcb := on_lhcb_normalized * 7000./nrj;')
