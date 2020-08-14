
from cc3d import CompuCellSetup
        
from JordanAmberModelSteppables import ODEModelSteppable
CompuCellSetup.register_steppable(steppable=ODEModelSteppable(frequency=1))

from JordanAmberModelSteppables import CellularModelSteppable
CompuCellSetup.register_steppable(steppable=CellularModelSteppable(frequency=1))

from JordanAmberModelSteppables import IFNPlotSteppable
CompuCellSetup.register_steppable(steppable=IFNPlotSteppable(frequency=1))

from JordanAmberModelSteppables import OutputSteppable
CompuCellSetup.register_steppable(steppable=OutputSteppable(frequency=1))

CompuCellSetup.run()
