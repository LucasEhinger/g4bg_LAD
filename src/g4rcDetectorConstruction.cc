#include "g4rcDetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "g4rcDetector.hh"
#include "g4rcMaterial.hh"
#include "G4SDManager.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Polycone.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4UserLimits.hh"
#include "G4SubtractionSolid.hh"

g4rcDetectorConstruction::g4rcDetectorConstruction() {

	fTarg = "H3";
	GetTargetIndex(fTarg);
	fHRS = "L";
	fHRSAngle = 17.5*deg;    

}


g4rcDetectorConstruction::~g4rcDetectorConstruction(){;}


G4VPhysicalVolume* g4rcDetectorConstruction::Construct() {
    

	fMaterial = g4rcMaterial::GetMaterialManager();

	G4double length = 25.0*cm;
	G4double radius = 0.65*cm;
	G4double thickness = 0.5*mm;

	G4Material* target_gas[4] = {fMaterial->H2_gas, fMaterial->D2_gas, fMaterial->T2_gas, fMaterial->He_gas};

	// WORLD LOGICAL
	
	double world_x, world_y, world_z;

	world_x = world_y = world_z = 1.5*m;

	G4Box* world_box = new G4Box("world",world_x,world_y,world_z);

	G4LogicalVolume* world_log
	= new G4LogicalVolume(world_box,fMaterial->air,"World",0,0,0);

	world_log->SetVisAttributes(G4VisAttributes::Invisible);

	// TARGET CHAMBER

	double t_chamber = 0.406*mm;
	double r_chamber = (1037./2.)*mm;
	
	G4Tubs* target_mother_tubs = new G4Tubs("targ_mother_tubs", 0., r_chamber+t_chamber, 25.*cm, 0.*deg, 360.*deg);
	G4LogicalVolume* target_mother_log = new G4LogicalVolume(target_mother_tubs,fMaterial->vacuum,"target_mother_logical",0,0,0);
	target_mother_log->SetVisAttributes(G4VisAttributes::Invisible);

	G4RotationMatrix* rotX_neg90 = new G4RotationMatrix();
	rotX_neg90->rotateX(-90.*deg);
	G4RotationMatrix* rotX_pos90 = new G4RotationMatrix();
	rotX_pos90->rotateX(90.*deg);

	// Target cell

	G4Tubs* cell_tubs = new G4Tubs("cell_tubs", 0., radius + thickness, length/2. + thickness, 0.*deg, 360.*deg);
	G4Tubs* cell_empt = new G4Tubs("cell_empt", 0., radius, length/2., 0.*deg, 360.*deg);

	G4SubtractionSolid* cell_sub = new G4SubtractionSolid("cell_sub", cell_tubs, cell_empt);
	G4LogicalVolume* cell_log = new G4LogicalVolume(cell_sub, fMaterial->aluminum, "cell_logical", 0, 0, 0);
	G4VPhysicalVolume* cell_phys = new G4PVPlacement(rotX_pos90,G4ThreeVector(), cell_log, "cell_physical", target_mother_log, false, 0);

	// Target gas volume
	
	G4Tubs* gas_tubs = new G4Tubs("gas_tubs", 0., radius, length/2., 0.*deg, 360.*deg);
//	G4LogicalVolume* gas_log = new G4LogicalVolume(gas_tubs, target_gas[fTargIndex], "gas_logical", 0,0,0);
	G4LogicalVolume* gas_log = new G4LogicalVolume(gas_tubs, fMaterial->H2_liquid, "gas_logical", 0,0,0);
	G4VisAttributes* gas_vis = new G4VisAttributes(G4Colour(0.,0.,1.));
	gas_log->SetVisAttributes(gas_vis);
	G4VPhysicalVolume* gas_phys = new G4PVPlacement(rotX_pos90,G4ThreeVector(), gas_log, "gas_physical", target_mother_log, false, 0);	


/*
	// Beryllium window (upstream of target)
	
	double t_be = 0.2003*mm;
	G4Tubs* be_tubs = new G4Tubs("be_tubs", 0., 5.*cm, t_be/2., 0.*deg, 360.*deg);
	G4LogicalVolume* be_log = new G4LogicalVolume(be_tubs, fMaterial->beryllium, "be_logical", 0,0,0);
	G4VisAttributes* be_vis = new G4VisAttributes(G4Colour(1.,1.,0.));
	be_log->SetVisAttributes(be_vis);
	G4VPhysicalVolume* be_phys = new G4PVPlacement(rotX_pos90,G4ThreeVector(0.,-30.*cm,0.), be_log, "be_physical", target_mother_log, false, 0);
*/
	
	G4VPhysicalVolume* target_mother_phys 
	= new G4PVPlacement(rotX_neg90,G4ThreeVector(), target_mother_log, "target_mother_physical", world_log, false, 0); 


	// Define the GEMs

	double w_gem1 = 10.*cm;
	double h_gem1 = 45.*cm;
	double w_gem2 = 13.*cm;
	double h_gem2 = 57.*cm;
	
	double t_gem = 0.5*cm;

	double r_gem1 = 70.*cm;
	double r_gem2 = 90.*cm;
	double gem_angle = 127.*deg;

	double x_gem1 = r_gem1*sin(gem_angle);
	double z_gem1 = r_gem1*cos(gem_angle);

	double x_gem2 = r_gem2*sin(gem_angle);
	double z_gem2 = r_gem2*cos(gem_angle);

	G4RotationMatrix* rot_gem = new G4RotationMatrix();
	rot_gem->rotateY(180.*deg - gem_angle);

	G4Box* gem1_box = new G4Box("gem1_box", w_gem1/2., h_gem1/2., t_gem/2.);  
	G4LogicalVolume* gem1_log = new G4LogicalVolume(gem1_box, fMaterial->vacuum, "gem1_log", 0, 0, 0);

	G4Box* gem2_box = new G4Box("gem2_box", w_gem2/2., h_gem2/2., t_gem/2.);  
	G4LogicalVolume* gem2_log = new G4LogicalVolume(gem2_box, fMaterial->vacuum, "gem2_log", 0, 0, 0);

	// Make the GEMs sensitive detectors

	G4SDManager* SDman = G4SDManager::GetSDMpointer();

	g4rcDetector* gem1_SD = new g4rcDetector("gem1_SD",101);
	SDman->AddNewDetector(gem1_SD);
	gem1_log->SetSensitiveDetector(gem1_SD);

	g4rcDetector* gem2_SD = new g4rcDetector("gem2_SD",102);
	SDman->AddNewDetector(gem2_SD);
	gem2_log->SetSensitiveDetector(gem2_SD);

	G4VPhysicalVolume* gem1_phys = new G4PVPlacement(rot_gem, G4ThreeVector(x_gem1, 0., z_gem1), gem1_log, "gem1_physical", world_log, false, 0); 
	G4VPhysicalVolume* gem2_phys = new G4PVPlacement(rot_gem, G4ThreeVector(x_gem2, 0., z_gem2), gem2_log, "gem2_physical", world_log, false, 0); 

	G4VPhysicalVolume* world_phys
	= new G4PVPlacement(0,G4ThreeVector(),world_log,"World",0,false,0);

	return world_phys;

}

void g4rcDetectorConstruction::GetTargetIndex(G4String targ) {

	fTargIndex = -1;
	G4String targList[4] = {"H1","H2","H3","He3"};

	for(int i = 0; i<4; i++) {
		if(targ == targList[i]) {
			fTargIndex = i;
		}		
	}

	if (fTargIndex < 0.) {
		G4cout << "Bad target selected!" << G4endl;
		exit(1);
	}
}
    
    
void g4rcDetectorConstruction::SetTarget(G4String targ) {
	fTarg = targ;
	GetTargetIndex(fTarg);
}






