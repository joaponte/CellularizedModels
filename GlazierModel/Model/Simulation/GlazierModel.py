from cc3d import CompuCellSetup

from GlazierModelSteppables import CellsInitializerSteppable

CompuCellSetup.register_steppable(steppable=CellsInitializerSteppable(frequency=1))

from GlazierModelSteppables import ImmuneRecruitmentSteppable

CompuCellSetup.register_steppable(steppable=ImmuneRecruitmentSteppable(frequency=1))

from GlazierModelSteppables import ViralSecretionSteppable

CompuCellSetup.register_steppable(steppable=ViralSecretionSteppable(frequency=1))

from GlazierModelSteppables import ImmuneCellKillingSteppable

CompuCellSetup.register_steppable(steppable=ImmuneCellKillingSteppable(frequency=1))

from GlazierModelSteppables import ChemotaxisSteppable

CompuCellSetup.register_steppable(steppable=ChemotaxisSteppable(frequency=1))

from GlazierModelSteppables import CytokineSecretionSteppable

CompuCellSetup.register_steppable(steppable=CytokineSecretionSteppable(frequency=1))

from GlazierModelSteppables import ViralLifeCycleSteppable

CompuCellSetup.register_steppable(steppable=ViralLifeCycleSteppable(frequency=1))

from GlazierModelSteppables import SimDataSteppable

CompuCellSetup.register_steppable(steppable=SimDataSteppable(frequency=1))

CompuCellSetup.run()
