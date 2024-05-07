#include "g4rcDetector.hh"
#include "G4SDManager.hh"
#include "G4String.hh"

g4rcDetector::g4rcDetector(G4String name, G4int detnum) : G4VSensitiveDetector(name)
{
	char colname[255];

	fDetNo = detnum;
	assert(fDetNo > 0);

	fHCID = -1;

	//    fTrackSecondaries = false;
	fTrackSecondaries = true;

	sprintf(colname, "genhit_%s_%d", name.data(), detnum);
	collectionName.insert(G4String(colname));
}

g4rcDetector::~g4rcDetector()
{
}

void g4rcDetector::Initialize(G4HCofThisEvent *)
{

	fHitColl = new g4rcDetectorHitsCollection(SensitiveDetectorName, collectionName[0]);
}

///////////////////////////////////////////////////////////////////////

G4bool g4rcDetector::ProcessHits(G4Step *step, G4TouchableHistory *)
{
	G4bool badedep = false;
	G4bool badhit = false;

	// Get touchable volume info
	G4TouchableHistory *hist =
		(G4TouchableHistory *)(step->GetPreStepPoint()->GetTouchable());
	G4int copyID = hist->GetReplicaNumber();

	G4StepPoint *prestep = step->GetPreStepPoint();
	G4Track *track = step->GetTrack();

	G4double Edep = step->GetTotalEnergyDeposit();
	G4int trackID = track->GetTrackID();
	G4int parentID = track->GetParentID();

	//    G4Material* material = track->GetMaterial();

	//    printf("Standard detector %d hit by %s!\n", fDetNo, track->GetParticleDefinition()->GetParticleName().data());

	g4rcDetectorHit *thishit = NULL;

	for (G4int i = fHitColl->entries() - 1; i >= 0; i--)
	{
		if (((*fHitColl)[i]->GetTrackID() == trackID))
		{
			// if ( ((*fHitColl)[i]->GetTrackID() == trackID) || ((*fHitColl)[i]->GetTrackID() == parentID)  ) {
			thishit = (*fHitColl)[i];
			break;
		}
	}

	if (track->GetDefinition()->GetPDGEncoding() != 2212) {
		//if the particle is not a proton, skip it
		badhit = true;
	}

	if (!thishit && !badhit)
	{
		thishit = new g4rcDetectorHit(fDetNo, copyID);
		thishit->SetTrackID(trackID);
		fHitColl->insert(thishit);
	}

	if (!badhit)
	{
		// Hit
		thishit->f3X = prestep->GetPosition();
		thishit->f3V = track->GetVertexPosition();
		thishit->f3D = track->GetVertexMomentumDirection();
		thishit->f3P = track->GetMomentum();

		thishit->fP = track->GetMomentum().mag();
		thishit->fE = track->GetTotalEnergy();
		thishit->fM = track->GetDefinition()->GetPDGMass();

		thishit->fTrID = track->GetTrackID();
		thishit->fmTrID = track->GetParentID();
		thishit->fPID = track->GetDefinition()->GetPDGEncoding();
		thishit->fTime = prestep->GetGlobalTime();

		if (Edep > 0.)
		{
			thishit->AddEdep(Edep);
		}

		// FIXME - Enumerate encodings
		thishit->fGen = (long int)track->GetCreatorProcess();
	}

	return !badedep && !badhit;
}

///////////////////////////////////////////////////////////////////////

void g4rcDetector::EndOfEvent(G4HCofThisEvent *HCE)
{
	G4SDManager *sdman = G4SDManager::GetSDMpointer();

	if (fHCID < 0)
	{
		fHCID = sdman->GetCollectionID(collectionName[0]);
	}

	HCE->AddHitsCollection(fHCID, fHitColl);
	return;
}
