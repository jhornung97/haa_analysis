import ROOT
from DataFormats.FWLite import Events, Handle
import numpy as np

events = Events('/ceph/bmaier/Haa/testsample/ZZTo2E2Nu_8B164E3B-0C2C-D74C-8335-2A05311FD80A.root')
handle = Handle('std::vector<pat::IsolatedTrack>')
label = ("isolatedTracks")

niso = 0
nel = 0
totnel = 0
pdgids = []

for event in events:
	event.getByLabel(label,handle)
	isos = handle.product()
	for iso in isos:
		if ((iso.pt() > 5 and (np.abs(iso.pdgId()) == 11 or np.abs(iso.pdgId()) == 13)) or iso.pt() > 10) and (np.abs(iso.pdgId()) < 15 or np.abs(iso.eta()) < 2.5) and ((np.abs(iso.dxy()) < 0.2 and np.abs(iso.dz()) < 0.1) or iso.pt() > 15) and ((iso.pfIsolationDR03().chargedHadronIso() < 5 and iso.pt() < 25) or iso.pfIsolationDR03().chargedHadronIso()/iso.pt() < 0.2):
			pdgids.append(np.abs(iso.pdgId()))
			if np.abs(iso.pdgId()) == 11:
				nel += 1
		if np.abs(iso.pdgId()) == 11:
			totnel += 1	
	niso += len(isos)

print('No. of isolated tracks: {}'.format(niso))
print('Total no. of electrons: {}'.format(totnel))
print('Isolated tracks that survived the cuts: {}'.format(len(pdgids)))
print('Electrons that survived: {}'.format(nel))
print('Unique PDG-Ids: {}'.format(np.unique(pdgids)))
