import os, socket
from kpf_cosmic_rays import __path__

RESULTSDIR = os.path.join(os.path.dirname(__path__[0]), 'results')
TABLEDIR = os.path.join(RESULTSDIR, 'tables')
DRIVERDIR = os.path.join(os.path.dirname(__path__[0]), 'drivers')

LOCALDIR = os.path.join(os.path.expanduser('~'), 'local', 'kpf_cosmic_rays')
if not os.path.exists(LOCALDIR):
    os.mkdir(LOCALDIR)

DATADIR = os.path.join(LOCALDIR, 'data')
