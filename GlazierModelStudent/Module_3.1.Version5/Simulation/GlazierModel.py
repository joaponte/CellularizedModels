from cc3d import CompuCellSetup

from GlazierModelSteppables import ModelSteppable

CompuCellSetup.register_steppable(steppable=ModelSteppable(frequency=1))

CompuCellSetup.run()
