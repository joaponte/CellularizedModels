
from cc3d import CompuCellSetup
        
from JordanAmberModelSteppables import ModelSteppable
CompuCellSetup.register_steppable(steppable=ModelSteppable(frequency=1))

from JordanAmberModelSteppables import PlaqueAssaySteppable
CompuCellSetup.register_steppable(steppable=PlaqueAssaySteppable(frequency=1))

CompuCellSetup.run()
