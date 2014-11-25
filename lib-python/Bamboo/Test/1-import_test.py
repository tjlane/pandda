#! .usr.localthonthon2.7.3-64bit.binthon

from Bamboo.Misc.SoakParser import *
from Bamboo.Misc.ModelHistory import *
from Bamboo.Analysis.Compare.Coords import *
from Bamboo.Analysis.Compare.Fragment import *
from Bamboo.Analysis.Compare.Angles import *
from Bamboo.Symmetry.Utils import *
from Bamboo.Maps.Utils import *
from Bamboo.Maps.Masking import *
from Bamboo.Maps.Convert import *
from Bamboo.Common.Command import *
from Bamboo.Common.Errors import *
from Bamboo.Common.File import *
from Bamboo.Common.Parser import *
from Bamboo.Common.Logs import *
from Bamboo.Common.Online.Utils import *
from Bamboo.Common.Path import *
from Bamboo.Utils.Ligand import *
from Bamboo.Utils.Constants import *
from Bamboo.Utils.PDB import *
from Bamboo.Rdkit.Smile import *
from Bamboo.Rdkit.Coords.Utils import *
from Bamboo.Rdkit.Mol import *
from Bamboo.Rdkit.Bonds.Dihedral import *
from Bamboo.Rdkit.Bonds.Fragment import *
from Bamboo.Utils.MTZ import *
from Bamboo.Utils.ProgramSelector import *
from Bamboo.Wrappers.PdbUtils import *
from Bamboo.Wrappers.Coot.Utils import *
from Bamboo.Wrappers.Coot.Scripts import *
from Bamboo.Wrappers.MtzUtils import *
from Bamboo.Wrappers.MapUtils import *
from Bamboo.Wrappers.Builders.Utils import *
from Bamboo.Wrappers.Builders.Grade import *
from Bamboo.Wrappers.Builders.Prodrg import *
from Bamboo.Wrappers.Builders.Writedict import *
from Bamboo.Wrappers.Builders.Elbow import *
from Bamboo.Wrappers.Refiners.Utils import *
from Bamboo.Wrappers.Refiners.Refmac import *
from Bamboo.Wrappers.BfactorUtils import *
from Bamboo.Wrappers.SymUtils import *
from Bamboo.Wrappers.Fitters.Utils import *
from Bamboo.Wrappers.Fitters.ARPwARP import *
from Bamboo.Wrappers.Fitters.Flynn import *
from Bamboo.Wrappers.Fitters.Rhofit import *
from Bamboo.Wrappers.Fitters.LigandFit import *
from Bamboo.Wrappers.CifUtils import *
from Bamboo.Macro.Utils import *
from Bamboo.Macro.Molecule import *
from Bamboo.Density.Edstats.Utils import *
from Bamboo.Density.Edstats.Score import *
from Bamboo.Xray.Utils import *
from Bamboo.Structure.Utils import *
from Bamboo.Structure.Contacts.Utils import *

print('SUCCESS: IMPORT ALL MODULES')

print ' ==> Test 1: {!s} failures.'.format(0)
