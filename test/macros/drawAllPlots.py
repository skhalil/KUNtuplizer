#!/usr/bin/env python

import subprocess, os
#from pdflatex import PDFLaTeX
#from tex import latex2pdf
options = [
 
    ['b', 'hCutflow_', True],
    ['notb', 'hCutflow_', True],
    ['b', 'hCutflow_eta_1p5_', True],
    ['notb', 'hCutflow_eta_1p5_', True],
    ['b', 'hCutflow_eta_2_', True],
    ['notb', 'hCutflow_eta_2_', True],
    ['b', 'hCutflow_eta_2p4_', True],
    ['notb', 'hCutflow_eta_2p4_', True],
    ['b', 'hCutflow_pt_5_', True],
    ['notb', 'hCutflow_pt_5_', True],
    ['b', 'hCutflow_pt_10_', True],
    ['notb', 'hCutflow_pt_10_', True],
    ['b', 'hCutflow_pt_20_', True],
    ['notb', 'hCutflow_pt_20_', True],

    ['b',    'hCutflow_pt_5_eta_1p5_', True], 
    ['notb', 'hCutflow_pt_5_eta_1p5_', True],
    ['b',    'hCutflow_pt_5_eta_2_',   True],
    ['notb', 'hCutflow_pt_5_eta_2_',   True], 
    ['b',    'hCutflow_pt_5_eta_2p4_', True],
    ['notb', 'hCutflow_pt_5_eta_2p4_', True],
    
    ['b',    'hCutflow_pt_10_eta_1p5_', True],
    ['notb', 'hCutflow_pt_10_eta_1p5_', True],
    ['b',    'hCutflow_pt_10_eta_2_',   True],
    ['notb', 'hCutflow_pt_10_eta_2_',   True], 
    ['b',    'hCutflow_pt_10_eta_2p4_', True],
    ['notb', 'hCutflow_pt_10_eta_2p4_', True],
        
    ['b',    'hCutflow_pt_20_eta_1p5_', True],
    ['notb', 'hCutflow_pt_20_eta_1p5_', True],
    ['b',    'hCutflow_pt_20_eta_2_',   True],
    ['notb', 'hCutflow_pt_20_eta_2_',   True], 
    ['b',    'hCutflow_pt_20_eta_2p4_', True],
    ['notb', 'hCutflow_pt_20_eta_2p4_', True],
    
    ['b', 'hSVPt_', False],
    ['notb', 'hSVPt_', False],

    ['b', 'hSVPt_f_', False],
    ['notb', 'hSVPt_f_', False],

    ['b', 'hSVPt_gen_', False],
    ['notb', 'hSVPt_gen_', False],

    ['b', 'hSVPt_gen_f_', False],
    ['notb', 'hSVPt_gen_f_', False],

    ['b', 'hSVEta_', False],
    ['notb', 'hSVEta_', False],

    ['b', 'hSVEta_f_', False],
    ['notb', 'hSVEta_f_', False],

    ['b', 'hSVEta_gen_', False],
    ['notb', 'hSVEta_gen_', False],

    ['b', 'hSVEta_gen_f_', False],
    ['notb', 'hSVEta_gen_f_', False],

    ['b', 'hSVMass_', False],
    ['notb', 'hSVMass_', False],

    ['b', 'hSVMass_f_', False],
    ['notb', 'hSVMass_f_', False],

    ['b', 'hSVMass_gen_', False],
    ['notb', 'hSVMass_gen_', False],

    ['b', 'hSVMass_gen_f_', False],
    ['notb', 'hSVMass_gen_f_', False]

    
    ]

command = 'python plot.py --isb={0:s} --hist={1:s} --cutflow' 

for option in options:

    s = command.format(
        option[0], option[1]
    )

    if option[2] == False:
        s = s.replace('--cutflow', '')
        print 'hello-->', s
    subprocess.call( ["echo --------------------------------------------------------------------------",""], shell=True)
    subprocess.call( ["echo %s"%s,""], shell=True)
    subprocess.call( ["echo --------------------------------------------------------------------------",""], shell=True)
    subprocess.call( [s, ""], shell=True )
    
    if option[2] == True:
        proc = subprocess.Popen(['pdflatex', 'Plots/'+option[1]+option[0]+'_table.tex'])
        proc.communicate()
        
        retcode = proc.returncode
        if not retcode == 0:
            os.unlink(option[1]+option[0]+'_table.pdf')
            raise ValueError('Error {} executing command: {}'.format(retcode, ' '.join(cmd)))

        os.rename(option[1]+option[0]+'_table.pdf', 'Plots/'+option[1]+option[0]+'_table.pdf')
        os.unlink('Plots/'+option[1]+option[0]+'_table.tex')
        os.unlink(option[1]+option[0]+'_table.log')
        os.unlink(option[1]+option[0]+'_table.aux')
