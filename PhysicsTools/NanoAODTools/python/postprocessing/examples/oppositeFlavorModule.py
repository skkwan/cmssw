# This is an example of a NanoAODTools Module to add one variable to nanoAODs.
# Note that:
# -the new variable will be available for use in the subsequent modules
# -it is possible to update the value for existing variables
#
# Example of using from command line:
# nano_postproc.py outDir /eos/cms/store/user/andrey/f.root -I PhysicsTools.NanoAODTools.postprocessing.examples.exampleModule exampleModuleConstr
#
# Example of running in a python script: see test/example_postproc.py
#

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from .helper import findLeadingPair, invariantMass, findLeadingPairFromTwoCollections
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from itertools import combinations

class exampleProducer(Module):
    def __init__(self, jetSelection, muoSelection, eleSelection):
        self.jetSel = jetSelection
        self.muoSel = muoSelection
        self.eleSel = eleSelection
        pass

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("m_muonele", "D")
        self.out.branch("nElectronsPassing", "I")
        self.out.branch("nMuonsPassing", "I")


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        # Initialize
        eventPasses = False 

        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        jets = Collection(event, "Jet")

        # First count the number of electrons, muons, and jets passing the baseline selections 
        filteredEle = list(filter(self.eleSel, electrons))
        filteredMuo = list(filter(self.muoSel, muons))
        filteredJet = list(filter(self.jetSel, jets))
        nElectronsPassing = len(filteredEle)
        nMuonsPassing = len(filteredMuo)
        nJetsPassing = len(filteredJet)

        # Initialize these branches
        m_muonele = 0

        print(f"nElectrons/Muons passing: {nElectronsPassing}, {nMuonsPassing}")

        # In the full processor we will do more advanced checks like cleaning the jets from the electrons and muons,
        # but at the minimum we need at least two jets. 
        if not ((nJetsPassing >= 2) and ((nElectronsPassing > 0) and (nMuonsPassing > 0))):
            return False 

        # Do combinations, if there are > 0 muons and >0 electrons in the event
        if ((nMuonsPassing > 0) and (nElectronsPassing > 0)):
            leadingOppFlavour = findLeadingPairFromTwoCollections(filteredMuo, filteredEle)
            m_muonele = invariantMass(leadingOppFlavour)
            print(f"Found leading opposite-flavour pair with mass {m_muonele}")

        # There must be >= 2 muons with OS and m_ll > 50, or >=2 electrons with OS and m_ee > 50
        if not ((nMuonsPassing > 0) and ((nElectronsPassing > 0) and (m_muonele > 50))):
            return False 

        # Fill branches 
        self.out.fillBranch("m_muonele", nMuonsPassing)
        self.out.fillBranch("nElectronsPassing", nElectronsPassing)
        self.out.fillBranch("nMuonsPassing", nMuonsPassing)
        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

oppositeFlavourModule = lambda: exampleProducer(jetSelection=lambda j: j.pt > 25, 
                                             # muons: see https://gitlab.cern.ch/mrherrma/zhmet/-/blob/main/src/chai/processors/tchizh.py#L327-338
                                             #            - TChiZH.py code says miniisotight, which is miniPFRelIso_all < 0.1
                                             muoSelection=lambda m: m.mediumId and (m.miniPFRelIso_all < 0.1) and (m.pt > 20) and (abs(m.eta) < 2.5) and (abs(m.ip3d) < 0.1) and (abs(m.dz) < 0.2),
                                             # electrons: see https://gitlab.cern.ch/mrherrma/zhmet/-/blob/main/src/chai/processors/tchizh.py#L343-362
                                             #            - Electron ID: Iso_WP90
                                             #            - Isolation: TChiZH.py code says iso of 0.1 to match miniisotight for muons, slides say 0.2
                                             #            - use 0.10 max ip3d (it's max 0.05 in the barrel and max 0.10 in the endcap)
                                             #            - use 0.2 max dz (it's max 0.1 in the barrel and max 0.2 in the endcap)
                                             eleSelection=lambda e: (e.mvaFall17V2Iso_WP90) and (e.miniPFRelIso_all < 0.1) and (e.pt > 20) and (abs(e.eta) < 2.5) and (abs(e.ip3d) < 0.10) and (abs(e.dz) < 0.2)
)

