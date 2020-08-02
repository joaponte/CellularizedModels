
from cc3d import CompuCellSetup
        

from GlazierModelSteppables import GlazierModelSteppable

CompuCellSetup.register_steppable(steppable=GlazierModelSteppable(frequency=1))


CompuCellSetup.run()
