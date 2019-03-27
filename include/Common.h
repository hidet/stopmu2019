#ifndef Common_h
#define Common_h 1

#include "G4SystemOfUnits.hh"
#include "globals.hh"
#include "TString.h"

static const G4String SAVEDIR=Form("./geant_data");

// --- acceptance mode ---                                               
//#define ACCEPTANCE        
//#define TARGETUNIFORM_GEANTINO
//#define TARGETUNIFORM

// ------ Number of Detectors ----
const G4int NbOfTES   = 256;
const G4int NbOfDEG   = 1;
const G4int NbOfTG    = 1;

// ---------- Max Multiplicity -----------
const G4int BeamMultiMax   = 1;
const G4int TESMultiMax    = 10;
const G4int TGMultiMax     = 10;
const G4int DEGMultiMax    = 10;

// -- Detector Construction --
// World (Vacuum)
const G4double world_size[3] = {10.0*m,10.0*m,10.0*m};
// Air
const G4double air_z    = 3.0*mm;
const G4double air_rmin = 0.0*mm;
const G4double air_rmax = 110.0*mm/2.;
const G4double air_pos[3]={0., 0., air_z/2.};
// accelerator window (Kapton)
const G4double aw_z    = 0.075*mm;//  75 um
const G4double aw_rmin = 0.0*mm;
const G4double aw_rmax = 110.0*mm/2.;
const G4double aw_pos[3]={0., 0., -0.1*mm};
// beam window (Kapton)
const G4double bw_z    = 0.100*mm;// 100 um
const G4double bw_rmin = 0.0*mm;
const G4double bw_rmax = 110.0*mm/2.;
const G4double bw_pos[3]={0., 0., air_pos[2]+air_z/2.+bw_z/2.};
// target (Ne gas)
const G4double target_z    = 702.8*mm;// 16.-3.+3.9+400+3.9+282
const G4double target_rmin = 0.0*mm;
const G4double target_rmax = 134.2*mm/2.;
const G4double target_pos[3]={0., 0., target_z/2.+bw_pos[2]+bw_z/2.};
// target side (SUS)
const G4double targetside1_z    = 159.75*mm;// 16.-3.+3.9+200.-114.3/2.
const G4double targetside1_rmin = 134.2*mm/2.;
const G4double targetside1_rmax = 139.8*mm/2.;
const G4double targetside1_pos[3]={0., 0., targetside1_z/2.+bw_pos[2]+bw_z/2.};
// target side (SUS)
const G4double targetside2_z    = 428.75*mm;//200.+3.9+282.-114.3/2.*mm
const G4double targetside2_rmin = 134.2*mm/2.;
const G4double targetside2_rmax = 139.8*mm/2.;
const G4double targetside2_pos[3]={0., 0., 16.+3.9+200.+114.3/2.+targetside2_z/2.};
// target side (SUS)
const G4double targetside3_z    = 114.3*mm;
const G4double targetside3_rmin = 134.2*mm/2.;
const G4double targetside3_rmax = 139.8*mm/2.;
const G4double targetside3_pos[3]={0., 0., 16.+3.9+200.*mm};
// target end (SUS)
const G4double targetend_z    = 2.0*mm;
const G4double targetend_rmin = 0.;
const G4double targetend_rmax = 139.8*mm/2.;
const G4double targetend_pos[3]={0., 0., target_pos[2]+target_z/2.+targetend_z/2.};
// 300K vacuum snout (SUS)
const G4double snout300k_rmin = 40.*mm/2.;
const G4double snout300k_rmax = 101.6*mm/2.;
const G4double snout300k_z    = 2.*mm;
const G4double snout300kwnd_rmin = 0.*mm;
const G4double snout300kwnd_rmax = 40.*mm/2.;
const G4double snout300kwnd_z    = 0.3*mm;// 300um
const G4double snout300kside_rmin = 97.4*mm/2.;
const G4double snout300kside_rmax = 101.6*mm/2.;
const G4double snout300kside_z    = 113.3*mm;
const G4double snout300k_pos[3]={70.*mm,0.,16.+3.9+200.*mm};
const G4double snout300kwnd_pos[3]={snout300k_pos[0],snout300k_pos[1],snout300k_pos[2]};
const G4double snout300kside_pos[3]={snout300k_pos[0]+snout300k_z/2.+snout300kside_z/2.,snout300k_pos[1],snout300k_pos[2]};
// 60K snout
const G4double snout60k_rmin = 30.*mm/2.;
const G4double snout60k_rmax = 85.*mm/2.;
const G4double snout60k_z    = 1.*mm;
const G4double snout60kwnd_rmin = 0.*mm;
const G4double snout60kwnd_rmax = 30.*mm/2.;
const G4double snout60kwnd_z    = 0.005*mm;// 5 um
const G4double snout60k_pos[3]={snout300k_pos[0]+3.,snout300k_pos[1],snout300k_pos[2]};
const G4double snout60kwnd_pos[3]={snout60k_pos[0],snout60k_pos[1],snout60k_pos[2]};
// 3K snout
const G4double snout3k_rmin = 25.*mm/2.;
const G4double snout3k_rmax = 70.*mm/2.;
const G4double snout3k_z    = 1.*mm;
const G4double snout3kwnd_rmin = 0.*mm;
const G4double snout3kwnd_rmax = 25.*mm/2.;
const G4double snout3kwnd_z    = 0.025*mm;// 25 um
const G4double snout3k_pos[3]={snout300k_pos[0]+6.,snout300k_pos[1],snout300k_pos[2]};
const G4double snout3kwnd_pos[3]={snout3k_pos[0],snout3k_pos[1],snout3k_pos[2]};
// 50mK snout
const G4double snout_rmin = 22.*mm/2.;
const G4double snout_rmax = 50.*mm/2.;
const G4double snout_z    = 1.*mm;
const G4double snoutwnd_rmin = 0.*mm;
const G4double snoutwnd_rmax = 22*mm/2.;
const G4double snoutwnd_z    = 0.025*mm;// 25 um
const G4double snout_pos[3]={snout300k_pos[0]+9.,snout300k_pos[1],snout300k_pos[2]};
const G4double snoutwnd_pos[3]={snout_pos[0],snout_pos[1],snout_pos[2]};
// TES Bi
const G4double calo_size[3] = {0.320*mm,0.300*mm,0.004*mm};
const G4double calo_pos_c[3]={snout300k_pos[0]+15.,snout300k_pos[1],snout300k_pos[2]};


#endif
