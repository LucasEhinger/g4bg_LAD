#include "LADDetectorConstructionGEMCreator.hh"
#include "g4rcMaterial.hh"
#include "g4rcDetector.hh"
#include "TString.h"
// #include "G4SBSDetectorConstruction.hh"
// #include "G4SBSGEMSD.hh"
// #include "sbstypes.hh"

#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UserLimits.hh"
#include "G4SDManager.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4TwoVector.hh"
#include "G4GenericTrap.hh"
#include "G4Polycone.hh"
#include "G4RotationMatrix.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

LADDetectorConstructionGEMCreator::LADDetectorConstructionGEMCreator()
{
    G4cout << "<LADDetectorConstructionGEMCreator> Info" << G4endl;
}

// This routine allows us to flexibly position GEM modules without code duplication:
void LADDetectorConstructionGEMCreator::BuildGEM(G4LogicalVolume *worldLV, G4SDManager *SDman)
//, G4RotationMatrix *rot, G4ThreeVector pos, unsigned int nplanes, vector<double> zplanes, vector<double> wplanes, vector<double> hplanes, G4String SDname )

{

    //   SDman = G4SDManager::GetSDMpointer();
    fMaterial = g4rcMaterial::GetMaterialManager();
    G4ThreeVector pos(0, 0, 0);
    G4int nplanes = 2;

    std::vector<double> zplanes;
    std::vector<double> wplanes;
    std::vector<double> hplanes;
    std::vector<int> layerid;
    G4double zspacing_SBStracker = 20.0 * cm;
    G4double zoffset_SBStracker = -20.0 * cm; // distance from the origin

    for (int i = 0; i < nplanes; i++)
    {
        zplanes.push_back(zoffset_SBStracker + i * zspacing_SBStracker);
        // wplanes.push_back( 60.0*cm );
        // hplanes.push_back( 200.0*cm );

        wplanes.push_back(120 * cm); // width
        hplanes.push_back(55 * cm);  // height
        layerid.push_back(101+i);//layerid
    }

    G4String SDname = "Harm/SBSGEM";

    // This routine will create and position a GEM tracker consisting of nplanes planes centered at position pos oriented with rotation rot wrt logical volume Mother.
    // The list of z coordinates, widths and heights of the planes are passed as arguments:

    G4bool GEMflip = FALSE;

    // How should we interpret the rotation Matrix rot? It is the rotation that orients the z axis of the tracker with respect to the mother volume.
    // Since pos is the nominal position of the tracker with respect to the mother volume,
    // the positioning of the centers of the tracker planes should be pos + zplane * tracker_zaxis
    G4ThreeVector zaxis(0, 0, 1);
    // zaxis *= rot->inverse();

    G4String TrackerPrefix = SDname;

    // Create sensitive detector for this tracker:
    G4String GEMSDname = SDname;
    G4String GEMSDname_nopath = SDname;
    GEMSDname_nopath.erase(0, SDname.rfind('/') + 1);
    G4String GEMcolname = GEMSDname_nopath; // We have to remove all the directory structure from the Hits collection name or else GEANT4 SDmanager routines will not handle correctly.
    GEMcolname += "HitsCollection";

    // This is related with the Sensitive Detector code. I need to consult how to work with it
    /*
      G4SBSGEMSD* GEMSD;

      if( !(GEMSD = (G4SBSGEMSD*) fDetCon->fSDman->FindSensitiveDetector(GEMSDname)) ){ //Make sure SD with this name doesn't already exist
      GEMSD = new G4SBSGEMSD( GEMSDname, GEMcolname );
      fDetCon->fSDman->AddNewDetector(GEMSD);
      (fDetCon->SDlist).insert(GEMSDname);
      fDetCon->SDtype[GEMSDname] = G4SBS::kGEM;
      }
    */

    //  fDetCon->InsertSDboundaryVolume( Mother->GetName(), GEMSDname );
    /*
      if( !( nplanes > 0 && zplanes.size() == nplanes && wplanes.size() == nplanes && hplanes.size() == nplanes ) ){
      G4cout << "MakeTracker called with invalid arguments, doing nothing..." << G4endl;
      return;
      }
    */
    // Best way to handle "flipping" of GEMs? Boxes are symmetric, just flip ordering of z positions!

    // Define z extent of various GEM layer components:

    const unsigned int nlayers = 24;

    double gempz[nlayers] = {
        120.0 * um, 3.0 * mm, 120.0 * um, 5.0 * um, 50.0 * um, 3.0 * mm,   // cover + cathode
        5.0 * um, 50.0 * um, 5.0 * um, 2.0 * mm,                           // GEM0
        5.0 * um, 50.0 * um, 5.0 * um, 2.0 * mm,                           // GEM1
        5.0 * um, 50.0 * um, 5.0 * um, 2.0 * mm,                           // GEM2
        10.0 * um, 50.0 * um, 180.0 * um, 120.0 * um, 3.0 * mm, 120.0 * um // readout + honeycomb
    };

    // Define GEM materials:

    G4Material *gempm[nlayers] = {
        fMaterial->NemaG10, fMaterial->NOMEX, fMaterial->NemaG10,
        fMaterial->Copper, fMaterial->Kapton, fMaterial->GEMgas,
        fMaterial->Copper, fMaterial->Kapton, fMaterial->Copper, fMaterial->GEMgas,
        fMaterial->Copper, fMaterial->Kapton, fMaterial->Copper, fMaterial->GEMgas,
        fMaterial->Copper, fMaterial->Kapton, fMaterial->Copper, fMaterial->GEMgas,
        fMaterial->Copper, fMaterial->Kapton, fMaterial->NemaG10, fMaterial->NemaG10, fMaterial->NOMEX, fMaterial->NemaG10};


    // G4Material *gempm[nlayers] = {
    //     fMaterial->Kapton, fMaterial->air, fMaterial->air,
    //     fMaterial->air, fMaterial->air, fMaterial->air,
    //     fMaterial->air, fMaterial->air, fMaterial->air, fMaterial->air,
    //     fMaterial->air, fMaterial->air, fMaterial->air, fMaterial->air,
    //     fMaterial->air, fMaterial->air, fMaterial->air, fMaterial->air,
    //     fMaterial->air, fMaterial->air, fMaterial->air, fMaterial->air, fMaterial->air, fMaterial->air};

    // G: NEMAG, B: NOMEX, Y: Cu, R: Kapton, W: GEMgas
    G4Colour GEMcolor[nlayers] = {G4Colour::Green(), G4Colour::Blue(), G4Colour::Green(),
                                  G4Colour::Yellow(), G4Colour::Red(), G4Colour::White(),
                                  G4Colour::Yellow(), G4Colour::Red(), G4Colour::Yellow(), G4Colour::White(),
                                  G4Colour::Yellow(), G4Colour::Red(), G4Colour::Yellow(), G4Colour::White(),
                                  G4Colour::Yellow(), G4Colour::Red(), G4Colour::Yellow(), G4Colour::White(),
                                  G4Colour::Yellow(), G4Colour::Red(), G4Colour::Green(), G4Colour::Green(), G4Colour::Blue(), G4Colour::Green()};

    unsigned int gidx, gpidx;

    // Determine the largest transverse dimensions of a plane in this tracker:
    double gemmaxw = 0.0;
    double gemmaxh = 0.0;

    for (gidx = 0; gidx < nplanes; gidx++)
    {
        if (wplanes[gidx] > gemmaxw)
            gemmaxw = wplanes[gidx];
        if (hplanes[gidx] > gemmaxh)
            gemmaxh = hplanes[gidx];
    }

    // Determine the total z extent of a plane in this tracker:
    double gempzsum = 0.0;
    for (gpidx = 0; gpidx < nlayers; gpidx++)
    {
        gempzsum += gempz[gpidx];
    }

    // implement aluminum shielding:
    // G4bool useAlshield = fDetCon->GetGEMuseAlshield();

    // G4Box *alshield_box;
    // G4LogicalVolume *alshield_log;

    G4double Al_thick = 0.0, Airgap_thick = 0.0;
    /*  if( useAlshield ){
    //G4cout << "building aluminum shielding for GEMs...";

    Al_thick = fDetCon->GetGEMAlShieldThick();
    Airgap_thick = fDetCon->GetGEMAirGapThick();


    //G4cout << " complete" << G4endl;
    }
    */

    char cgidx[255];
    char cgpidx[255];

    G4Box *gpbox;
    G4LogicalVolume *gplog;

    G4VisAttributes *gemvisatt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));

    gemvisatt->SetForceWireframe(true);

    G4VisAttributes *gemsdvisatt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));

    for (gidx = 0; gidx < nplanes; gidx++)
    {
        sprintf(cgidx, "%02d", gidx);

        G4String gemboxname = TrackerPrefix + G4String("_gembox_");
        gemboxname += cgidx;
        G4String gemlogname = TrackerPrefix + G4String("_gemlog_");
        gemlogname += cgidx;

        G4Box *gembox = new G4Box(gemboxname, wplanes[gidx] / 2.0, hplanes[gidx] / 2.0, gempzsum / 2.0);
        // G4LogicalVolume *gemlog = new G4LogicalVolume(gembox, fMaterial->air, gemlogname, 0, 0, 0); // orig
        G4LogicalVolume *gemlog = new G4LogicalVolume(gembox, fMaterial->vacuum, gemlogname, 0, 0, 0);

        gemlog->SetVisAttributes(gemvisatt);

        double ztemp = 0.0;
        double sign = 1.0;
        if (GEMflip)
        { // this reverses ordering of GEM materials front-to-back!
            ztemp = gempzsum;
            sign = -1.0;
        }
        /*
        if( useAlshield ){ //put Aluminum shielding upstream of GEM:
          G4String shieldboxname = TrackerPrefix + G4String("alshield_gem")   + cgidx + G4String("box");
          G4String shieldlogname = TrackerPrefix + G4String("alshield_gem") + cgidx + G4String("log");

          G4String frontshieldphysname = TrackerPrefix + G4String("alshieldfront_gem") + cgidx + G4String("phys");
          G4String backshieldphysname = TrackerPrefix + G4String("alshieldback_gem") + cgidx + G4String("phys");

          G4cout << "Building and placing GEM aluminum shielding, plane = " << gidx << G4endl;

          G4Box *alshield_box = new G4Box( shieldboxname, wplanes[gidx]/2.0, hplanes[gidx]/2.0, Al_thick/2.0 );
          G4LogicalVolume *alshield_log = new G4LogicalVolume( alshield_box,  Materials->Aluminum, shieldlogname, 0, 0, 0 );

          G4VisAttributes *alshield_visatt = new G4VisAttributes( G4Colour( 0.6, 0.6, 0.6 ) );
          //alshield_visatt->SetForceWireframe( true );

          alshield_log->SetVisAttributes( alshield_visatt );

          //position the front shielding upstream of the GEM with possible air gap:
          G4ThreeVector frontshieldpos = pos + (zplanes[gidx] - gempzsum/2.0 - Airgap_thick - Al_thick/2.0) * zaxis;

          new G4PVPlacement( 0, frontshieldpos, alshield_log, frontshieldphysname, Mother, true, gidx+1, false );

          //position the back shielding immediately downstream of the GEM with no air gap:
          G4ThreeVector backshieldpos = frontshieldpos + (gempzsum + Al_thick + Airgap_thick ) * zaxis;

          new G4PVPlacement( 0, backshieldpos, alshield_log, backshieldphysname, Mother, true, gidx+1, false );
        }
        */
        G4double angleW = (90 - 30 - 23) * deg;

        G4RotationMatrix rGEM[2];
        rGEM[0].rotateY(-(23 + 30) * deg);
        rGEM[1].rotateY(-(23 + 30) * deg);

        G4double distanceW = 70 * cm;

        G4double vSeparationX = distanceW * cos(angleW);
        G4double vSeparationY = 0;
        G4double vSeparationZ = -distanceW * sin(angleW);

        G4double wallSep = zspacing_SBStracker;

        G4double vSpaceX = wallSep * cos(angleW);
        G4double vSpaceY = 0;
        G4double vSpaceZ = -wallSep * sin(angleW);

        for (gpidx = 0; gpidx < nlayers; gpidx++)
        {
            sprintf(cgpidx, "_%02d_%03d_", gidx, gpidx);

            G4String gempboxname = TrackerPrefix + G4String("_gemplane") + cgpidx + G4String("box");
            G4String gemplogname = TrackerPrefix + G4String("_gemplane") + cgpidx + G4String("log");
            G4String gempphysname = TrackerPrefix + G4String("_gemplane") + cgpidx + G4String("phy");

            ztemp += sign * gempz[gpidx] / 2.0;

            gpbox = new G4Box(gempboxname, wplanes[gidx] / 2.0, hplanes[gidx] / 2.0, gempz[gpidx] / 2.0);
            // gplog = new G4LogicalVolume(gpbox, gempm[gpidx], gemplogname, 0, 0, 0); // orig
            gplog = new G4LogicalVolume(gpbox, fMaterial->vacuum, gemplogname, 0, 0, 0);

            new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, ztemp - gempzsum / 2.0), gplog, gempphysname, gemlog, false, 0, false);

            ztemp += sign * gempz[gpidx] / 2.0;

            // Assign sensitive volume: why 5?  // SPR: This is the gas drift region
            if (gpidx == 5)
            {
                g4rcDetector *GEM_SD = new g4rcDetector(Form("GEM_SD_%i", layerid[gidx]), layerid[gidx]);
                SDman->AddNewDetector(GEM_SD);
                gplog->SetSensitiveDetector(GEM_SD);
                // gplog->SetSensitiveDetector(GEMSD);

                gplog->SetVisAttributes(G4VisAttributes(GEMcolor[gpidx])); // every logical volume has a color (useless but funny)

                // Until we implement actual strips/wires in the GEM construction, the detmap is irrelevant for the GEMs
            }
            else
            {
                gplog->SetVisAttributes(G4VisAttributes::GetInvisible());
            }
        }

        G4String gemname = TrackerPrefix + G4String("_gemphys_") + cgidx;
        // Now place the fully constructed GEM plane in the mother logical volume:
        // G4ThreeVector plane_pos = pos + G4ThreeVector( 0.0, 0.0, zplanes[gidx] );
        //    G4ThreeVector plane_pos = pos + zplanes[gidx] * zaxis;
        // Now we are positioning the GEM AFTER positioning all of its components inside its logical volume:

        G4ThreeVector plane_pos = G4ThreeVector(vSeparationX, vSeparationY, vSeparationZ) +
                                  G4ThreeVector(vSpaceX * gidx, vSpaceY * gidx, vSpaceZ * gidx);

        //  The position vector should be modified and not following the tracker
        // I will do an approach, similar to the Hodo

        new G4PVPlacement(G4Transform3D(rGEM[gidx], plane_pos), gemlog, gemname, worldLV, true, gidx + 1, false);
    }
}

LADDetectorConstructionGEMCreator::~LADDetectorConstructionGEMCreator()
{
    G4cout << "<~LADDetectorConstructionGEMCreator> --> DONE" << G4endl;
}
