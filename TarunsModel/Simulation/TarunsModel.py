
from cc3d import CompuCellSetup

from TarunsModelSteppables import TarunsModelSteppable
CompuCellSetup.register_steppable(steppable=TarunsModelSteppable(frequency=1))

from TarunsModelSteppables import ChemotaxisSteppable
CompuCellSetup.register_steppable(steppable=ChemotaxisSteppable(frequency=1))

from TarunsModelSteppables import PlotsSteppable
CompuCellSetup.register_steppable(steppable=PlotsSteppable(frequency=1))

CompuCellSetup.run()
