
from cc3d import CompuCellSetup
        
from AmberFluModelSteppables import AmberFluModelSteppable
CompuCellSetup.register_steppable(steppable=AmberFluModelSteppable(frequency=1))

from CoinfectionModelSteppables import CellularModelSteppable
CompuCellSetup.register_steppable(steppable=CellularModelSteppable(frequency=1))

from CoinfectionModelSteppables import Data_OutputSteppable
CompuCellSetup.register_steppable(steppable=Data_OutputSteppable(frequency=1))

from CoinfectionModelSteppables import PlotSteppable
CompuCellSetup.register_steppable(steppable=PlotSteppable(frequency=1))

CompuCellSetup.run()
