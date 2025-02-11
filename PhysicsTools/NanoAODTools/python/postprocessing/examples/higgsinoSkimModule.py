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
from .helper import findLeadingPair, invariantMass
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
        self.out.branch("nMuonsPassing", "I")
        self.out.branch("nElectronsPassing", "I")
        self.out.branch("m_dimuon", "D")
        self.out.branch("m_diele", "D")

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
        m_dimuon = 0
        m_diele = 0

        # In the full processor we will do more advanced checks like cleaning the jets from the electrons and muons,
        # but at the minimum we need at least two jets. 
        if not ((nJetsPassing >= 2) and ((nElectronsPassing >= 2) or (nMuonsPassing >= 2))):
            return False 

        # Do combinations, if there are two or more muons
        if (nMuonsPassing >= 2):
            leadingMuons = findLeadingPair(filteredMuo)
            m_dimuon = invariantMass(leadingMuons)
            print(f"Found leading muon pair with pT {leadingMuons[0].pt} and {leadingMuons[1].pt}, total mass of {m_dimuon}")

        # Do combinations if there are two or more electrons
        if (nElectronsPassing >= 2):
            leadingElectrons = findLeadingPair(filteredEle)
            m_diele = invariantMass(leadingElectrons)
            print(f"Found leading electron pair with pT {leadingElectrons[0].pt} and {leadingElectrons[1].pt}, total mass of {m_diele}")

        # There must be >= 2 muons with OS and m_ll > 50, or >=2 electrons with OS and m_ee > 50
        if not (((nMuonsPassing >= 2) and (m_dimuon > 50)) or ((nElectronsPassing >= 2) and (m_diele > 50))):
            return False 

        # Fill branches 
        self.out.fillBranch("nMuonsPassing", nMuonsPassing)
        self.out.fillBranch("m_dimuon", m_dimuon)
        self.out.fillBranch("nElectronsPassing", nElectronsPassing)
        self.out.fillBranch("m_diele", m_diele)

        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

higgsinoSkimModule = lambda: exampleProducer(jetSelection=lambda j: j.pt > 25, 
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

