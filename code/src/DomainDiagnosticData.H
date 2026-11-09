#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef DOMAIN_DIAGNOSTIC_DATA_H
#define DOMAIN_DIAGNOSTIC_DATA_H

#include "CH_HDF5.H"
#include "Vector.H"
#include "LevelSigmaCS.H"
#include "ParmParse.H"

#include "NamespaceHeader.H"

/// Storage and I/O for some common domain wide diagnostics against time
/** 
    Currently intended for inclusion in the CF-friendly hdf5 I/O.

    Includes 'Scalars' from Table A3 from Nowicki et al, Geosci. Model Dev., 9, 4521–4545, 2016 
    ---------------------------------------------------------------------------------
  
    Tier 2
    ------
    land_ice_mass (kg) [VOLUME, m^3]
    land_ice_mass_not_displacing_sea_water (kg) [VAF, m^3]
    
    Tier 3
    ------
    grounded_land_ice_area (m^2)
    floating_ice_shelf_area (m^2)
    tendency_of_land_ice_mass_due_to_surface_mass_balance (kg s^-1) [SMB, m/a]
    tendency_of_land_ice_mass_due_to_basal_mass_balance (kg s^-1) [BMB, m/a]
    tendency_of_land_ice_mass_due_to_calving (kg s^-1) [m/a]

*/


class AmrIce;

namespace MaskedIntegration
{
  void integrateScalarInside
  (Real& a_integral, 
   function<bool(Real h, Real f, int mask)> a_inside,
   const Vector<RefCountedPtr<LevelSigmaCS > >& a_coords,
   const Vector<LevelData<FArrayBox>* >& a_integrand, 
   const Vector<LevelData<FArrayBox>* >& a_topography,
   const Vector<LevelData<FArrayBox>* >& a_thickness,
   const Vector<LevelData<FArrayBox>* >& a_iceFrac,
   const Vector<LevelData<FArrayBox>* >& a_sectorMaskFraction,
   const Vector<Real>& a_dx, const Vector<int>& a_ratio,
   int a_finestLevel, 
   int a_maskNo, int a_maskComp);

  void integrateDischargeInside
  (Real& a_sumDischarge, Real& a_sumDivUH,
   function<bool(Real h, Real f, int mask)> a_inside,
   const Vector<RefCountedPtr<LevelSigmaCS > >& a_coords,
   const Vector<LevelData<FluxBox>* >& a_fluxOfIce,
   const Vector<LevelData<FArrayBox>* >& a_topography,
   const Vector<LevelData<FArrayBox>* >& a_thickness,
   const Vector<LevelData<FArrayBox>* >& a_iceFrac,
   const Vector<Real>& a_dx, const Vector<int>& a_ratio, 
   const Vector<LevelData<FArrayBox>* >& a_sectorMaskFraction,
   int a_finestLevel,  
   int a_maskNo, int a_maskComp);
  
}


#define HDF5_SUBGROUP_NAME "bisicles_domain_diagnostic_data"
#define CFIO_DIAGNOSTIC_TIME_NAME "time"

#define CFIO_DIAGNOSTIC_ICE_VOLUME_SHORT_NAME "lim"
#define CFIO_DIAGNOSTIC_ICE_VAF_SHORT_NAME "limnsw"
#define CFIO_DIAGNOSTIC_GROUNDED_ICE_AREA_SHORT_NAME "iareagr"
#define CFIO_DIAGNOSTIC_FLOATING_ICE_AREA_SHORT_NAME "iareafl"
#define CFIO_DIAGNOSTIC_UPPER_SURFACE_FLUX_SHORT_NAME "tendacabf"
#define CFIO_DIAGNOSTIC_LOWER_SURFACE_FLUX_SHORT_NAME "tendlibmassbf"
#define CFIO_DIAGNOSTIC_LOWER_SURFACE_GROUNDED_ICE_FLUX_SHORT_NAME "tendlibmassbfgr"
#define CFIO_DIAGNOSTIC_LOWER_SURFACE_FLOATING_ICE_FLUX_SHORT_NAME "tendlibmassbffl"
#define CFIO_DIAGNOSTIC_CALVING_FLUX_SHORT_NAME "tendlicalvf"
#define CFIO_DIAGNOSTIC_ICE_FRONT_CALVING_AND_MELTING_SHORT_NAME "tendlifmassbf"
#define CFIO_DIAGNOSTIC_GL_FLUX_SHORT_NAME "tendligroundf"

#define CFIO_DIAGNOSTIC_ICE_VOLUME_CF_NAME "land_ice_mass"
#define CFIO_DIAGNOSTIC_ICE_VAF_CF_NAME "land_ice_mass_not_displacing_sea_water"
#define CFIO_DIAGNOSTIC_GROUNDED_ICE_AREA_CF_NAME "grounded_ice_sheet_area"
#define CFIO_DIAGNOSTIC_FLOATING_ICE_AREA_CF_NAME "floating_ice_sheet_area"
#define CFIO_DIAGNOSTIC_UPPER_SURFACE_FLUX_CF_NAME "tendency_of_land_ice_mass_due_to_surface_mass_balance"
#define CFIO_DIAGNOSTIC_LOWER_SURFACE_FLUX_CF_NAME "tendency_of_land_ice_mass_due_to_basal_mass_balance"
#define CFIO_DIAGNOSTIC_LOWER_SURFACE_GROUNDED_ICE_FLUX_CF_NAME "tendency_of_land_ice_mass_due_to_basal_mass_balance_floating"
#define CFIO_DIAGNOSTIC_LOWER_SURFACE_FLOATING_ICE_FLUX_CF_NAME "tendency_of_land_ice_mass_due_to_basal_mass_balance_grounded"
#define CFIO_DIAGNOSTIC_CALVING_FLUX_CF_NAME "tendency_of_land_ice_mass_due_to_calving"
#define CFIO_DIAGNOSTIC_GL_FLUX_CF_NAME "tendency_of_grounded_ice_mass_due_to_discharge"
#define CFIO_DIAGNOSTIC_ICE_FRONT_CALVING_AND_MELTING_CF_NAME "tendency_of_land_ice_mass_due_to_calving_and_ice_front_melting"

#define CFIO_DIAGNOSTIC_ICE_VOLUME_LONG_NAME "Ice mass"
#define CFIO_DIAGNOSTIC_ICE_VAF_LONG_NAME "Ice mass not displaced by seawater"
#define CFIO_DIAGNOSTIC_GROUNDED_ICE_AREA_LONG_NAME "Area covered by grounded ice"
#define CFIO_DIAGNOSTIC_FLOATING_ICE_AREA_LONG_NAME "Area covered by floating ice"
#define CFIO_DIAGNOSTIC_UPPER_SURFACE_FLUX_LONG_NAME "Total SMB flux"
#define CFIO_DIAGNOSTIC_LOWER_SURFACE_FLUX_LONG_NAME "Total BMB flux"
#define CFIO_DIAGNOSTIC_LOWER_SURFACE_GROUNDED_ICE_FLUX_LONG_NAME "Total BMB flux beneath grounded ice"
#define CFIO_DIAGNOSTIC_LOWER_SURFACE_FLOATING_ICE_FLUX_LONG_NAME "Total BMB flux beneath floating ice"
#define CFIO_DIAGNOSTIC_CALVING_FLUX_LONG_NAME "Total calving flux"
#define CFIO_DIAGNOSTIC_GL_FLUX_LONG_NAME "Total grounding line flux"
#define CFIO_DIAGNOSTIC_ICE_FRONT_CALVING_AND_MELTING_LONG_NAME "Total calving and ice front melting flux"

struct cfDiagnostic
{
  std::string short_name;
  std::string long_name;
  std::string cf_name;
  std::string units;
  Vector<Real>* data;
};

class DomainDiagnosticData
{
  Vector<Real> m_time; // time 
  Vector<Real> m_ice_vaf; // total volume above flotation
  Vector<Real> m_ice_volume; // total ice volume
  Vector<Real> m_ice_grounded_area; // grounded ice area
  Vector<Real> m_ice_floating_area; // floating ice area
  Vector<Real> m_ice_total_smb; // total SMB over ice covered regions
  Vector<Real> m_ice_total_bmb; // total BMB over ice covered regions
  Vector<Real> m_ice_grounded_total_bmb; // total BMB over floating ice
  Vector<Real> m_ice_floating_total_bmb; // total BMB over floating ice
  Vector<Real> m_ice_total_calving_flux; // total calving flux
  Vector<Real> m_ice_total_gl_flux; // total gl flux
  Vector<Real> m_ice_total_calving_and_ice_front_melting_flux; // total calving and ice front melting flux
  Vector<Real> m_error_ice_mass_con;

  Real m_previous_time, m_previous_volume; // used for computing the calving rate
  
  void setCFdata();
  
public:
  /// default constructor
  DomainDiagnosticData();

private:
  // disable copy constructor
  DomainDiagnosticData(const DomainDiagnosticData&);
  // disable assignment 
  DomainDiagnosticData& operator=(const DomainDiagnosticData&);
public:

  void initDiagnostics(AmrIce& a_amrIce,
		       const Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
		       const Vector<DisjointBoxLayout>& a_grids,
		       const Vector<int>& a_refRatio, Real a_crseDx,
		       Real a_time, int a_finestLevel);

  /// record diagnostics given the current state of a_amrIce
  void record(AmrIce& a_amrIce);
  /// reset 
  void reset(int a_len=0);
  
  // integrate all quantities over time, return a DiagnosticData containing the means
  DomainDiagnosticData* timeMean();
  
#ifdef CH_USE_HDF5
  /// write to hdf5
  void write(HDF5Handle& a_handle);
  /// read from hdf5
  void read(HDF5Handle& a_handle);
#endif

  /// store domain diagnostics to output later
  /**
     [0] time, [1] volume of ice, [2] volume above flotation, [3] area of grounded ice,
     [4] area of floating ice, [5] total SMB, [6] total BMB, 
     [7] total BMB beneath floating ice, [8] total calving flux,
     [9] total calving and ice front melting flux.  

   */
  //Real m_diagnostic_values[16]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, .0}; 
  //cfDiagnostic m_cf_stuff[10]; 
  vector<cfDiagnostic> m_cf_stuff;
  
  /// if true, report sum(grounded ice) as well as sum(total ice)
  bool m_report_grounded_ice;

  /// if true, report area of grounded ice and area of floating ice
  bool m_report_area;

  /// if true, report sum(basal flux) and sum(surface flux)
  bool m_report_total_flux;

  /// if true, report sum(melt from calving model) and sum(calved ice)
  bool m_report_calving;

  /// diagnostic quantity -- sum(ice) at initial time
  Real m_initialSumIce;

  /// diagnostic quantity -- sum(ice) at last timestepa
  Real m_lastSumIce;

  /// diagnostic quantity -- sum(calved ice) at last timestep
  Real m_lastSumCalvedIce;

  /// diagnostic quantity -- sum(grounded ice) at initial time
  Real m_initialSumGroundedIce;

  /// diagnostic quantity -- sum(grounded ice) at last timestep
  Real m_lastSumGroundedIce;

  /// diagnostic quantity -- sum(ice over flotation) at initial time
  Real m_initialVolumeAboveFlotation;
 
  /// diagnostic quantity -- sum(ice over flotation) at last timestepa
  Real m_lastVolumeAboveFlotation;

  /// if true, apply a mask for the domain diagnostic data output (to be read in)
  bool m_mask_domain_diagnostic_data;

  /// for the domain diagnostic data mask
  std::string m_domain_diagnostic_data_mask_file;
  int m_num_domain_diagnostic_data_masks;
  Vector<std::string> m_domain_diagnostic_data_mask_names;
  RefCountedPtr<LevelData<FArrayBox> > m_domain_diagnostic_data_mask;
  Real m_domain_diagnostic_data_mask_dx;
};

  
#include "NamespaceFooter.H"

#endif


