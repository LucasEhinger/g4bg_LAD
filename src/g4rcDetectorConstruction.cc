#include "g4rcDetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "g4rcDetector.hh"
#include "g4rcScintDetector.hh"
#include "g4rcMaterial.hh"
#include "G4SDManager.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Box.hh"
#include "G4Trap.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Polycone.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4OpticalSurface.hh"
#include "G4SubtractionSolid.hh"
#include "G4VisAttributes.hh"

g4rcDetectorConstruction::g4rcDetectorConstruction() {

	fTarg = "T2";
	fHRS = "L";

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

	double world_x, world_y, world_z;
	double det_x, det_y, det_z;

	world_x = world_y = world_z = 275*cm;
	det_x = 15*cm;
	det_y = 6*cm;
	det_z = 6*cm;

	G4Box* world_box = new G4Box("World",world_x,world_y,world_z);

	G4LogicalVolume* world_log
	= new G4LogicalVolume(world_box,fMaterial->air,"World",0,0,0);

	world_log->SetVisAttributes(G4VisAttributes::Invisible);

	G4VPhysicalVolume* world_phys
	= new G4PVPlacement(0,G4ThreeVector(),world_log,"World",0,false,0);


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
	G4VPhysicalVolume* cell_phys = new G4PVPlacement(0,G4ThreeVector(), cell_log, "cell_physical", world_log, false, 0); 

	// Target gas volume
	
	G4Material* target_gas;

	if(fTarg == "T2") {
		target_gas = fMaterial->T2_gas;
	}

	G4Tubs* gas_tubs = new G4Tubs("gas_tubs", 0., radius, length/2., 0.*deg, 360.*deg);
	G4LogicalVolume* gas_log = new G4LogicalVolume(gas_tubs, target_gas, "gas_logical", 0,0,0);
	G4VisAttributes* gas_vis = new G4VisAttributes(G4Colour(0.,0.,1.));
	gas_log->SetVisAttributes(gas_vis);
	G4VPhysicalVolume* gas_phys = new G4PVPlacement(0,G4ThreeVector(), gas_log, "gas_physical", world_log, false, 0);	

	// Beryllium window (upstream of target)
	
	double t_be = 0.2003*mm;
	G4Tubs* be_tubs = new G4Tubs("be_tubs", 0., 5.*cm, t_be/2., 0.*deg, 360.*deg);
	G4LogicalVolume* be_log = new G4LogicalVolume(be_tubs, fMaterial->beryllium, "be_logical", 0,0,0);
	G4VisAttributes* be_vis = new G4VisAttributes(G4Colour(1.,1.,0.));
	be_log->SetVisAttributes(be_vis);
	G4VPhysicalVolume* be_phys = new G4PVPlacement(0,G4ThreeVector(0.,0.,-50.*cm), be_log, "be_physical", world_log, false, 0);
	

	// Rotation matrix needed for next two volumes
	G4RotationMatrix* rotx90 = new G4RotationMatrix();
	rotx90->rotateX(-90.*deg);

	// Target chamber exit window
	double t_ch = 0.406*mm;
	double r_ch = (1037./2.)*mm;
	G4Tubs* chamber_tubs = new G4Tubs("chamber_tubs", r_ch, r_ch+t_ch, 15.*cm ,0.*deg, 180.*deg);
	G4LogicalVolume* chamber_log = new G4LogicalVolume(chamber_tubs, fMaterial->aluminum, "chamber_logical", 0,0,0);
	G4VPhysicalVolume* chamber_phys = new G4PVPlacement(rotx90,G4ThreeVector(), chamber_log, "chamber_physical", world_log, false, 0); 

	// Q1 entrance window
	double t_q1w = 0.305*mm;
	double d_ch_q1 = 81.6*cm;
	G4Tubs* q1w_tubs = new G4Tubs("q1w_tubs", r_ch+t_ch+d_ch_q1, r_ch+t_ch+d_ch_q1+t_q1w, 25.*cm ,0.*deg, 180.*deg);
	G4LogicalVolume* q1w_log = new G4LogicalVolume(q1w_tubs, fMaterial->kapton, "q1w_logical", 0,0,0);
	G4VPhysicalVolume* q1w_phys = new G4PVPlacement(rotx90,G4ThreeVector(), q1w_log, "q1w_physical", world_log, false, 0); 




	return world_phys;

}


    
    
    
    


