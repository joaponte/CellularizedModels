
from cc3d import CompuCellSetup

from JordanIFNModelSteppables import ODEModelSteppable
CompuCellSetup.register_steppable(steppable=ODEModelSteppable(frequency=1))

from JordanIFNModelSteppables import CellularModelSteppable
CompuCellSetup.register_steppable(steppable=CellularModelSteppable(frequency=1))

from JordanIFNModelSteppables import PlotODEModelSteppable
CompuCellSetup.register_steppable(steppable=PlotODEModelSteppable(frequency=1))

CompuCellSetup.run()
