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

g4rcDetectorConstruction::g4rcDetectorConstruction() {

	fTarg = "T2";
	fHRS = "L";
	fHRSAngle = 17.5*deg;

	if(fHRS=="L") {
		fHRSAngle*=-1.0;
	}

	length = 25.0*cm;
	radius = 0.65*cm;

	if(fTarg == "T2") {
	entrance_window = 0.253*mm;
	exit_window = 0.343*mm;
		if(fHRS == "L") {
			entrance_wall = 0.473*mm;
			mid_wall = 0.435*mm;
			exit_wall = 0.379*mm;
		} else if(fHRS == "R") {
			entrance_wall = 0.425*mm;
			mid_wall = 0.447*mm;
			exit_wall = 0.406*mm;
		}
	}
		
    
}


g4rcDetectorConstruction::~g4rcDetectorConstruction(){;}


G4VPhysicalVolume* g4rcDetectorConstruction::Construct() {
    

	fMaterial = g4rcMaterial::GetMaterialManager();



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

	// Target chamber exit window
	G4Tubs* chamber_tubs = new G4Tubs("chamber_tubs", r_chamber, r_chamber+t_chamber, 20.*cm ,0.*deg, 360.*deg);
	G4LogicalVolume* chamber_log = new G4LogicalVolume(chamber_tubs, fMaterial->aluminum, "chamber_logical", 0,0,0);
	G4VPhysicalVolume* chamber_phys = new G4PVPlacement(0,G4ThreeVector(), chamber_log, "chamber_physical", target_mother_log, false, 0); 

	// Target cell
	const int nPoints = 10;
	double cell_z[nPoints], cell_r[nPoints];
	
	// outer edge
	cell_z[0] = -((length/2.) + entrance_window) ;
	cell_r[0] = 0.;

	cell_z[1] = -((length/2.) + entrance_window);
	cell_r[1] = radius + entrance_wall;

	cell_z[2] = 0.;
	cell_r[2] = radius + mid_wall;

	cell_z[3] = (length/2.) + exit_window;
	cell_r[3] = radius + exit_wall;   
 
	cell_z[4] = (length/2.) + exit_window;
	cell_r[4] = 0.;

	//inner edge
	cell_z[5] = length/2.;
	cell_r[5] = 0.;

	cell_z[6] = length/2.;
	cell_r[6] = radius;
	
	cell_z[7] = 0.;
	cell_r[7] = radius; 

	cell_z[8] = -length/2.;
	cell_r[8] = radius;  

	cell_z[9] = -length/2.;
	cell_r[9] = 0.;
	

	G4Polycone* cell_polycone = new G4Polycone("cell_polycone", 0.*deg, 360.*deg, nPoints, cell_r, cell_z);
	G4LogicalVolume* cell_log = new G4LogicalVolume(cell_polycone, fMaterial->aluminum, "cell_logical",0,0,0);
	G4VPhysicalVolume* cell_phys = new G4PVPlacement(rotX_pos90,G4ThreeVector(), cell_log, "cell_physical", target_mother_log, false, 0); 

	// Target gas volume
	
	G4Material* target_gas;

	if(fTarg == "T2") {
		target_gas = fMaterial->T2_gas;
	}

	G4Tubs* gas_tubs = new G4Tubs("gas_tubs", 0., radius, length/2., 0.*deg, 360.*deg);
	G4LogicalVolume* gas_log = new G4LogicalVolume(gas_tubs, target_gas, "gas_logical", 0,0,0);
	G4VisAttributes* gas_vis = new G4VisAttributes(G4Colour(0.,0.,1.));
	gas_log->SetVisAttributes(gas_vis);
	G4VPhysicalVolume* gas_phys = new G4PVPlacement(rotX_pos90,G4ThreeVector(), gas_log, "gas_physical", target_mother_log, false, 0);	

	double maxStep = 5.*cm;
	fStepLimit = new G4UserLimits(maxStep);
	gas_log->SetUserLimits(fStepLimit);

	// Beryllium window (upstream of target)
	
	double t_be = 0.2003*mm;
	G4Tubs* be_tubs = new G4Tubs("be_tubs", 0., 5.*cm, t_be/2., 0.*deg, 360.*deg);
	G4LogicalVolume* be_log = new G4LogicalVolume(be_tubs, fMaterial->beryllium, "be_logical", 0,0,0);
	G4VisAttributes* be_vis = new G4VisAttributes(G4Colour(1.,1.,0.));
	be_log->SetVisAttributes(be_vis);
	G4VPhysicalVolume* be_phys = new G4PVPlacement(rotX_pos90,G4ThreeVector(0.,-30.*cm,0.), be_log, "be_physical", target_mother_log, false, 0);
	

	G4VPhysicalVolume* target_mother_phys 
	= new G4PVPlacement(rotX_neg90,G4ThreeVector(), target_mother_log, "target_mother_physical", world_log, false, 0); 


	// Q1 entrance window
	double t_q1_window = 0.305*mm;
	double r_q1 = 15.*cm;
	double d_q1 = r_chamber + t_chamber + 81.6*cm;

	double z_q1 = d_q1*cos(fHRSAngle);
	double x_q1 = -d_q1*sin(fHRSAngle);

	G4RotationMatrix* rot_HRS = new G4RotationMatrix();
	rot_HRS->rotateY(fHRSAngle);
	
	G4Tubs* q1_tubs = new G4Tubs("q1_tubs", 0., r_q1, t_q1_window/2., 0.*deg, 360.*deg);
	G4LogicalVolume* q1_log = new G4LogicalVolume(q1_tubs, fMaterial->kapton, "q1_log", 0,0,0);
	G4VPhysicalVolume* q1_phys = new G4PVPlacement(rot_HRS,G4ThreeVector(x_q1,0.,z_q1), q1_log, "q1_physical", world_log, false, 0);  
	
	G4VPhysicalVolume* world_phys
	= new G4PVPlacement(0,G4ThreeVector(),world_log,"World",0,false,0);

	// Detector
	
	G4Tubs* det_tubs = new G4Tubs("det_tubs",0.,20.*cm,1.*mm,0.*deg,360.*deg);
	G4LogicalVolume* det_log = new G4LogicalVolume(det_tubs, fMaterial->vacuum, "det_log", 0,0,0);
		
	G4SDManager* SDman = G4SDManager::GetSDMpointer();
	g4rcDetector* detSD = new g4rcDetector("detSD",101);
	SDman->AddNewDetector(detSD);
	det_log->SetSensitiveDetector(detSD);

	double z_det = (d_q1+5.*cm)*cos(fHRSAngle);
	double x_det = -(d_q1+5.*cm)*sin(fHRSAngle);

	G4VPhysicalVolume* det_phys = new G4PVPlacement(rot_HRS, G4ThreeVector(x_det,0.,z_det), det_log, "det_physical", world_log, false, 0);

	return world_phys;

}


    
    
    
    


