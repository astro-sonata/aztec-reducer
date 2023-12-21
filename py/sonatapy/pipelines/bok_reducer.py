'''
A command line wrapper for bok testing
'''
import os, glob
from warnings import warn
from sonatapy import Bok
from sonatapy.exceptions import MissingDataException
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from multiprocessing import Pool
import numpy as np

def reorganize_datadir(datadir, standard):
    '''
    Properly reorganize the data directory for running the pipeline code
    '''
    
    if os.path.exists(os.path.join(datadir, 'flats')):
        # Assume that all of the subdirectories exist then
        warn('Not rearranging directory! I hope you formatted it correctly!')
        return # we don't need to reorganize

    # move the obvious ones
    subdirs = ['flats', 'focus', 'zero', 'reduced', standard] 
    for subdir in subdirs:        
        outdir = os.path.join(datadir, subdir)
        os.mkdir(outdir)

        if subdir == 'flats':
            file_ext = 'cont'
        else:
            file_ext = subdir
        
        files = glob.glob(os.path.join(datadir, '*'+file_ext+'*'))
        for f in files:
            outpath = os.path.join(outdir, os.path.basename(f))
            if os.path.basename(f) != file_ext:
                print(f'Moving {f} to {outpath}')
                os.rename(f, outpath)
    
    # now assume the ones that are left are all targets
    targ_files = glob.glob(os.path.join(datadir, '*.fits'))
    targ_names = np.unique([f.split('.')[0].split('_')[0] for f in targ_files])
    for targ in targ_names:
        outdir = os.path.join(datadir, targ)
        os.mkdir(outdir)
        files = glob.glob(os.path.join(datadir, '*'+os.path.basename(targ)+'*'))
        for f in files:
            outpath = os.path.join(outdir, os.path.basename(f))
            if f != targ:
                print(f'Moving {f} to {outpath}')
                os.rename(f, outpath)
    
            
def reduce1(kwargs):
    '''
    Just reduce one of the objects
    '''
    print()
    print(kwargs)
    print()

    try:
        bok_data = Bok(**kwargs)
    except MissingDataException as e:
        # we can just skip this one
        warn(f'Skipping! {e}')
        return

    outfile, wave, spec1d = bok_data.reduce_data(debug=False, overwrite=True)

def main():

    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('--datadir', '-d', required=True, help="The name of the directory with the data")
    p.add_argument('--standard_name', '-s', required=True, help="The name of the standard star")
    p.add_argument('--ncores', '-n', default=os.cpu_count(), help="The number of cores to run on")
    args = p.parse_args()

    # prep the data directory
    reorganize_datadir(args.datadir, args.standard_name)

    # now we can run the pipeline
    subdirs = set(glob.glob(os.path.join(args.datadir, '*')))
    rm = {os.path.join(args.datadir, args.standard_name),
          os.path.join(args.datadir, 'flats'),
          os.path.join(args.datadir, 'focus'),
          os.path.join(args.datadir, 'zero'),
          os.path.join(args.datadir, 'reduced')
          }

    obj_names = list(subdirs - rm)

    # package the arguments
    red_in = [{'obj_name': os.path.basename(o),
              'standard_name': args.standard_name,
              'datadir': args.datadir}
              for o in obj_names
              ]

    n = int(args.ncores)
    if n == 1 or n == 0:
        for kwarg in red_in:
            reduce1(kwarg)
    else:
        with Pool(n) as p:
            p.map(reduce1, red_in)
    
    
if __name__ == '__main__':
    main()
