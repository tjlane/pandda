import os, sys, glob

#######################################

if __name__=='__main__':

    if len(sys.argv) > 1:
        proc_dirs = sorted([os.path.abspath(p) for p in sys.argv[1:] if os.path.exists(p)])
    else:
        proc_dirs = sorted([os.path.abspath(p) for p in glob.glob('./*') if os.path.isdir(p)])

    print 'Processing {} Directories'.format(len(proc_dirs))

    for d in proc_dirs:

        print 'Processing {}'.format(d)
        os.chdir(d)
        if os.path.exists('./refine.pdb'):
            continue
        else:
            os.system('giant.quick_refine params=phenix.params program=phenix *.cif *-ensemble-model.pdb *-pandda-input.mtz')

