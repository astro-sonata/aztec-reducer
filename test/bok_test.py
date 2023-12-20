'''
A command line wrapper for bok testing
'''
import os
from sonatapy import Bok
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():

    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('--obj')
    args = p.parse_args()
    
    obj_name = args.obj #'sn2023tpl' #'at2019qiz'
    standard_name = 'feige110'
    datadir = os.path.join(os.path.dirname(os.path.dirname(os.getcwd())), 'bok', '2023Dec14')

    bok_data = Bok(datadir, obj_name, standard_name)
    wave, spec1d = bok_data.reduce_data(debug=True)
    
    fig, ax = plt.subplots(figsize=(18,6))
    ax.plot(wave, spec1d)
    ax.set_ylabel(f'Flux [{spec1d.unit}]')
    ax.set_xlabel(r'Wavelength [$\AA$]')
    fig.savefig('reduced_spectrum.png')
    
if __name__ == '__main__':
    main()
