
from cc3d import CompuCellSetup
        
from SaenzIFNModelSteppables import SaenzModelSteppable
CompuCellSetup.register_steppable(steppable=SaenzModelSteppable(frequency=1))

from SaenzIFNModelSteppables import CellularModelSteppable
CompuCellSetup.register_steppable(steppable=CellularModelSteppable(frequency=1))

from SaenzIFNModelSteppables import StatisticsSteppable
CompuCellSetup.register_steppable(steppable=StatisticsSteppable(frequency=1))

CompuCellSetup.run()
