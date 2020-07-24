
from cc3d import CompuCellSetup
        
from JordanAmberModelSteppables import ModelSteppable
CompuCellSetup.register_steppable(steppable=ModelSteppable(frequency=1))

from JordanAmberModelSteppables import IFNPlotSteppable
CompuCellSetup.register_steppable(steppable=IFNPlotSteppable(frequency=1))

CompuCellSetup.run()
