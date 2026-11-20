#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "SurfaceFlux.H"
#include "ComplexSurfaceFlux.H"
#include "LevelDataSurfaceFlux.H"
#include "ISMIP6OceanForcing.H"
#include "GroundingLineLocalizedFlux.H"
#include "HotspotFlux.H"
#include "BuelerGIA.H"
#include <map>
#ifdef HAVE_PYTHON
#include "PythonInterface.H"
#endif
#include "IceConstants.H"
#include "FineInterp.H"
#include "CoarseAverage.H"
#include "computeNorm.H"
#include "BisiclesF_F.H"
#include "ParmParse.H"
#include "AmrIceBase.H"
#include "FortranInterfaceIBC.H"
#include "FillFromReference.H"
#include "ReadLevelData.H"

#include "NamespaceHeader.H"

  /// factory method
  /** return a pointerto a new SurfaceFlux object
   */

SurfaceFlux* 
zeroFlux::new_surfaceFlux()
{
  zeroFlux* newPtr = new zeroFlux;
  return static_cast<SurfaceFlux*>(newPtr);
}

  /// define source term for thickness evolution and place it in flux
  /** dt is included in case one needs integrals or averages over a
      timestep
  */
void
zeroFlux::surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
			       const AmrIceBase& a_amrIce, 
			       int a_level, Real a_dt)
{
  DataIterator dit = a_flux.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      a_flux[dit].setVal(0.0);
    }
}


constantFlux::constantFlux() : m_isValSet(false)
{
}

SurfaceFlux* 
constantFlux::new_surfaceFlux()
{
  constantFlux* newPtr = new constantFlux;
  newPtr->m_fluxVal = m_fluxVal;
  newPtr->m_isValSet = m_isValSet;
  return static_cast<SurfaceFlux*>(newPtr);
}

  /// define source term for thickness evolution and place it in flux
  /** dt is included in case one needs integrals or averages over a
      timestep
  */
void
constantFlux::surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				   const AmrIceBase& a_amrIce, 
				   int a_level, Real a_dt)
{
  CH_assert(m_isValSet);
  DataIterator dit = a_flux.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      a_flux[dit].setVal(m_fluxVal);
    }
}


///
void
constantFlux::setFluxVal(const Real& a_fluxVal) 
{
  m_fluxVal = a_fluxVal; 
  // input value is in meters/year divide by secondsperyear 
  // to get flux in meters/second
  // slc: switch to flux in m/a
  //m_fluxVal/= secondsperyear;
  
  m_isValSet = true;
}


SurfaceFlux* SurfaceFlux::parse(const char* a_prefix)
{
  
  SurfaceFlux* ptr = NULL;
  std::string type = "";
  
  ParmParse pp(a_prefix);
  pp.query("type",type);
  
  if (type == "zeroFlux")
    {
      ptr = new zeroFlux;
    }
  else if (type == "constantFlux")
    {
      constantFlux* constFluxPtr = new constantFlux;
      Real fluxVal;
      pp.get("flux_value", fluxVal);
      constFluxPtr->setFluxVal(fluxVal);
      ptr = static_cast<SurfaceFlux*>(constFluxPtr);
    }
  else if (type == "hotspotFlux")
    {
      HotspotFlux* hotspotFluxPtr = new HotspotFlux;
      Real fluxVal;
      pp.get("flux_value", fluxVal);
      hotspotFluxPtr->setFluxVal(fluxVal);
      Vector<Real> vect(SpaceDim,0.0);

      pp.getarr("radius",vect,0,SpaceDim);
      RealVect radius(D_DECL(vect[0], vect[1],vect[2]));      

      pp.getarr("center",vect,0,SpaceDim);
      RealVect center(D_DECL(vect[0], vect[1],vect[2]));

      hotspotFluxPtr->setSpotLoc(radius, center);
      
      Real startTime = -1.2345e300;
      Real stopTime = 1.2345e300;
      pp.query("start_time", startTime);
      pp.query("stop_time", stopTime);
      hotspotFluxPtr->setSpotTimes(startTime, stopTime);
      
      ptr = static_cast<SurfaceFlux*>(hotspotFluxPtr);
    }
  else if (type == "LevelData")
    {
      std::string fileFormat;
      pp.get("fileFormat",fileFormat);
      int n;
      pp.get("n",n);
      int offset = 0;
      pp.query("offset",offset);

      Real startTime = 0.0, timeStep = 1.0;
      pp.query("startTime", startTime);
      pp.query("timeStep", timeStep);
      std::string name = "flux";
      pp.query("name", name);
      bool linearInterp = true;
      pp.query("linearInterp", linearInterp);

      Real defaultValue = 0.0;
      pp.query("defaultValue", defaultValue);
      
      RefCountedPtr<std::map<Real,std::string> > tf
	(new std::map<Real,std::string>);
      
      for (int i =0; i < n; i++)
	{
	  char* file = new char[fileFormat.length()+32];
	  sprintf(file, fileFormat.c_str(),i + offset);
	  tf->insert(make_pair(startTime + Real(i)*timeStep, file));
	  delete[] file;
	}
      
      LevelDataSurfaceFlux* ldptr = new LevelDataSurfaceFlux(tf,name,linearInterp, defaultValue);
      ptr = static_cast<SurfaceFlux*>(ldptr);
    }
  else if (type == "MultiLevelData")
    {
      std::string fileFormat;
      pp.get("fileFormat",fileFormat);
      int n;
      pp.get("n",n);
      int offset = 0;
      pp.query("offset",offset);

      Real startTime = 0.0, timeStep = 1.0;
      pp.query("startTime", startTime);
      pp.query("timeStep", timeStep);
      std::string name = "flux";
      pp.query("name", name);
      bool linearInterp = true;
      pp.query("linearInterp", linearInterp);

      Real defaultValue = 0.0;
      pp.query("defaultValue", defaultValue);
      
      RefCountedPtr<std::map<Real,std::string> > tf
	(new std::map<Real,std::string>);
      
      for (int i =0; i < n; i++)
	{
	  char* file = new char[fileFormat.length()+32];
	  sprintf(file, fileFormat.c_str(),i + offset);
	  tf->insert(make_pair(startTime + Real(i)*timeStep, file));
	  delete[] file;
	}
      
      MultiLevelDataSurfaceFlux* ldptr = new MultiLevelDataSurfaceFlux(tf,name,linearInterp, defaultValue);
      ptr = static_cast<SurfaceFlux*>(ldptr);
    }
  else if (type == "fortran")
    {
      // don't have the context here to actually set values, but
      // we can at least allocate the object here and return it 
      fortranInterfaceFlux* fifptr = new fortranInterfaceFlux;
      ptr = static_cast<SurfaceFlux*>(fifptr);
    }
  else if (type == "piecewiseLinearFlux")
    {
      int n = 1;  
      pp.query("n",n);
      Vector<Real> vabs(n,0.0);
      Vector<Real> vord(n,0.0);
      pp.getarr("abscissae",vabs,0,n);
      pp.getarr("ordinates",vord,0,n);
      
      Real dmin = -1.0;
      pp.query("minWaterDepth",dmin);
      
      PiecewiseLinearFlux* pptr = new PiecewiseLinearFlux(vabs,vord,dmin);
      ptr = static_cast<SurfaceFlux*>(pptr);
    }
  else if (type == "depthPowerFlux")
    {
      int nComponents = 1;
	  pp.query("nComponents",nComponents);
      Vector<Real> power(nComponents);  
	  pp.queryarr("power",power,0,power.size());
      Vector<Real> coef(nComponents);
	  pp.queryarr("coefficient",coef,0,coef.size());
      Real dmin = -1.0;
      pp.query("minWaterDepth",dmin);
      Real dtc = 0.0;
      pp.query("thermoclineDepth",dtc);
      bool multiplySqrtWaterDepth = false;
      pp.query("multiplySqrtWaterDepth",multiplySqrtWaterDepth);
      
      DepthPowerFlux* pptr = new DepthPowerFlux(power,coef,dmin,dtc,multiplySqrtWaterDepth,nComponents);
      ptr = static_cast<SurfaceFlux*>(pptr);
    }
  else if (type == "productFlux")
    {
      std::string flux1Prefix(a_prefix);
      flux1Prefix += ".flux1";
      SurfaceFlux* flux1Ptr = parse(flux1Prefix.c_str());
      if (flux1Ptr == NULL)
	{
      pout() << flux1Prefix << "has undefined flux1" << std::endl;
	  MayDay::Error("undefined flux1 in productFlux");
	}

      std::string flux2Prefix(a_prefix);
      flux2Prefix += ".flux2";
      SurfaceFlux* flux2Ptr = parse(flux2Prefix.c_str());
      if (flux2Ptr == NULL)
	{
	  MayDay::Error("undefined flux2 in productFlux");
	}
 

      ptr = static_cast<SurfaceFlux*>
	(new ProductSurfaceFlux(flux1Ptr->new_surfaceFlux(),
				flux2Ptr->new_surfaceFlux()));

      
      delete flux1Ptr;
      delete flux2Ptr;
    }
  else if (type == "conservedSumIceShelfFlux")
  {
      std::string prefix(a_prefix);
      prefix += ".flux";
      SurfaceFlux* fptr = parse(prefix.c_str());
      if (fptr == NULL)
	{
	  fptr = new zeroFlux;
	}
      
      ParmParse ppg("geometry");
      int iter(0);
      ppg.query("compute_ocean_connection_iter", iter);
      bool check_ocean_connected(false);
      if (iter > 0) pp.query("check_ocean_connected", check_ocean_connected);
      bool conserve_sum(true);
      pp.query("conserve_sum",conserve_sum);
      Real known_sum(0.0);
      bool sum_known = pp.contains("known_sum");
      if (sum_known) pp.get("known_sum", known_sum);

      ptr = static_cast<SurfaceFlux*>(new ConservedSumIceShelfFlux(fptr, check_ocean_connected, 
			      conserve_sum, sum_known, known_sum));
 
  }
  else if (type == "maskedFlux")
    {
      std::string groundedPrefix(a_prefix);
      groundedPrefix += ".grounded";
      SurfaceFlux* groundedPtr = parse(groundedPrefix.c_str());
      if (groundedPtr == NULL)
	{
	  groundedPtr = new zeroFlux;
	}
 
      std::string floatingPrefix(a_prefix);
      floatingPrefix += ".floating";
      SurfaceFlux* floatingPtr = parse(floatingPrefix.c_str());
      if (floatingPtr == NULL)
	{
	  floatingPtr = new zeroFlux;
	}

      std::string openLandPrefix(a_prefix);
      openLandPrefix += ".openLand";
      SurfaceFlux* openLandPtr = parse(openLandPrefix.c_str());
      if (openLandPtr == NULL)
	{
	  openLandPtr = groundedPtr->new_surfaceFlux();
	}

      
      std::string openSeaPrefix(a_prefix);
      openSeaPrefix += ".openSea";
      SurfaceFlux* openSeaPtr = parse(openSeaPrefix.c_str());
      if (openSeaPtr == NULL)
	{
	  openSeaPtr = floatingPtr->new_surfaceFlux();
	}


      bool floating_check_ocean_connected = false;
      pp.query("floating_check_ocean_connected",floating_check_ocean_connected);
      
      ptr = static_cast<SurfaceFlux*>
	(new MaskedFlux(groundedPtr->new_surfaceFlux(),
			floatingPtr->new_surfaceFlux(),
			openSeaPtr->new_surfaceFlux(),
			openLandPtr->new_surfaceFlux(),
			floating_check_ocean_connected) );
      
      delete groundedPtr;
      delete floatingPtr;
      delete openSeaPtr;
      delete openLandPtr;
    }
  else if (type == "boxBoundedFlux")
    {
      Vector<Real> tmp(SpaceDim); 
      pp.getarr("lo",tmp,0,SpaceDim);
      RealVect lo (D_DECL(tmp[0], tmp[1],tmp[2]));
      pp.getarr("hi",tmp,0,SpaceDim);
      RealVect hi (D_DECL(tmp[0], tmp[1],tmp[2]));

      Vector<Real> time(2);
      time[0] = -1.2345678e+300;
      time[1] = 1.2345678e+300;
      pp.queryarr("time",time,0,2);
           
      std::string prefix(a_prefix);
      prefix += ".flux";
      SurfaceFlux* fluxPtr = parse(prefix.c_str());
      CH_assert(fluxPtr != NULL);
      BoxBoundedFlux bbf(lo,hi,time[0],time[1],fluxPtr);
      ptr = static_cast<SurfaceFlux*>(bbf.new_surfaceFlux());

    }
  else if (type == "axbyFlux")
   {
     Real a; 
     pp.get("a",a);
     
     std::string xpre(a_prefix);
     xpre += ".x";
     SurfaceFlux* x = parse(xpre.c_str());
     
     Real b; 
     pp.get("b",b);
     
     std::string ypre(a_prefix);
     ypre += ".y";
     SurfaceFlux* y = parse(ypre.c_str());
    
     AxbyFlux axbyFlux(a,x,b,y);
     ptr = static_cast<SurfaceFlux*>(axbyFlux.new_surfaceFlux());

   }
  else if (type == "compositeFlux")
   {
     
     
     int nElements;
     pp.get("nElements",nElements);
     
     std::string elementPrefix(a_prefix);
     elementPrefix += ".element";

     Vector<SurfaceFlux*> elements(nElements);
     for (int i = 0; i < nElements; i++)
       {
	 std::string prefix(elementPrefix);
	 char s[32];
	 sprintf(s,"%i",i);
	 prefix += s;
	 ParmParse pe(prefix.c_str());
	 elements[i] = parse(prefix.c_str());
	 CH_assert(elements[i] != NULL);
	 
       }
     CompositeFlux compositeFlux(elements);
     ptr = static_cast<SurfaceFlux*>(compositeFlux.new_surfaceFlux());
   
   }
  else if (type == "normalizedFlux")
   {
     
     Real amplitude;
     pp.get("amplitude",amplitude);
     std::string prefix(a_prefix);
     prefix += ".direction";
     SurfaceFlux* direction = parse(prefix.c_str());
     if (direction == NULL)
       {
	 pout() << " no flux defined in " << prefix  << std::endl;
	 MayDay::Error("no flux defined in <normalized flux>.direction");
       }
     NormalizedFlux flux(direction, amplitude);
     ptr = static_cast<SurfaceFlux*>(flux.new_surfaceFlux());
   
   }
  else if (type == "targetThicknessFlux")
   {
     
     Real timescale = 1.0e+4;
     pp.get("timescale",timescale);
     std::string prefix(a_prefix);
     prefix += ".target";
     SurfaceFlux* target = parse(prefix.c_str());
     if (target == NULL)
       {
	 pout() << " no flux defined in " << prefix  << std::endl;
	 MayDay::Error("no flux defined in <target flux>.target");
	 
       }
     
     TargetThicknessFlux flux(target, timescale);
     ptr = static_cast<SurfaceFlux*>(flux.new_surfaceFlux());
   
   }
  else if (type == "groundingLineLocalizedFlux")
    {
      Real powerOfThickness = 0.0;
      pp.query("powerOfThickness",powerOfThickness);

      std::string glPrefix(a_prefix);
      glPrefix += ".groundingLine";
      SurfaceFlux* glPtr = parse(glPrefix.c_str());
      if (glPtr == NULL)
	{
	  glPtr = new zeroFlux;
	}
       
      std::string ambientPrefix(a_prefix);
      ambientPrefix += ".ambient";
      SurfaceFlux* ambientPtr = parse(ambientPrefix.c_str());
      if (ambientPtr == NULL)
	{
	  ambientPtr = new zeroFlux;
	}
      
      ptr = static_cast<SurfaceFlux*>
	(new GroundingLineLocalizedFlux(glPtr->new_surfaceFlux(),
					ambientPtr->new_surfaceFlux(),
					powerOfThickness ));
	 
      delete glPtr;
      delete ambientPtr;

    }
  else if (type == "FloatingDivUHLocalizedFlux")
    {
      std::string prefix(a_prefix);
      prefix +=  ".flux";
      SurfaceFlux* fluxptr = parse(prefix.c_str());

      //ParmParse pp2(prefix.c_str());
      Real dx;
      pp.get("mesh_spacing",dx);
      
      ptr = static_cast<SurfaceFlux*>
	(new  FloatingDivUHLocalizedFlux(fluxptr->new_surfaceFlux(), dx));
      delete fluxptr;
      
    }
  else if (type == "IMSIP6OceanForcing")
    {
      ptr = new ISMIP6OceanForcing(pp);
    }
#if FFTW_3
  else if (type == "buelerGIA") {
    // Read and set material constants.
    ParmParse ppCon("constants");
    // Defaults copied from AmrIce.cpp
    Real  m_iceDensity = 910.0; 
    Real  m_seaWaterDensity = 1028.0;
    Real  m_mantleDensity = 3313.0; 
    Real  m_gravity = 9.81;
    ppCon.query("ice_density",m_iceDensity);
    ppCon.query("gravity",m_gravity);
    ppCon.query("mantle_density",m_mantleDensity);
    ppCon.query("seaWater_density",m_seaWaterDensity);
    BuelerGIAFlux* buelerFlux = new BuelerGIAFlux(m_iceDensity, m_mantleDensity, m_gravity, m_seaWaterDensity);

    // Domain and FFT properties.
    ParmParse ppAmr("amr");
    Vector<int> ancells(3);
    ppAmr.getarr("num_cells",ancells, 0, ancells.size());
    int Nx, Ny;
    Nx = ancells[0];
    Ny = ancells[1];

    // see if the domain is offset from the origin
    Vector<int> offset(SpaceDim,0);
    IntVect domainOffset;
    ppAmr.queryarr("domainLoIndex", offset, 0, SpaceDim);
    domainOffset[0] = offset[0];
    domainOffset[1] = offset[1];

    ParmParse ppMain("main");
    Vector<Real> domsize(3);
    ppMain.getarr("domain_size", domsize, 0, 3);
    Real Lx, Ly;
    Lx = domsize[0];
    Ly = domsize[1];
    int pad;
    pp.get("pad", pad);
    bool m_inside_box = true;
    pp.query("inside_box", m_inside_box);
    int m_gia_box_lox = -9999999;
    pp.query("box_lox", m_gia_box_lox);
    int m_gia_box_hix = 9999999;
    pp.query("box_hix", m_gia_box_hix);
    int m_gia_box_loy = -9999999;
    pp.query("box_loy", m_gia_box_loy);
    int m_gia_box_hiy = 9999999;
    pp.query("box_hiy", m_gia_box_hiy);
    buelerFlux->setDomain(Nx, Ny, Lx, Ly, domainOffset, pad,
                          m_inside_box,
                          m_gia_box_lox, m_gia_box_hix,
                          m_gia_box_loy, m_gia_box_hiy); 

    Real t = 0.;
    ppAmr.query("offsetTime", t);
    buelerFlux->setUpdatedTime(t);

    int nlayers;
    pp.get("nlayers", nlayers);

    // Interpret the number of layers and type the properties accordingly.
    if ( nlayers == 1) {
      Real visc;
      pp.get("visc", visc);
      buelerFlux->setViscosity(visc, 1);
    }
    else if ( nlayers == 2 ) {
      Vector<Real> visc(2);
      pp.getarr("visc",visc,0,nlayers);
      Real thk;
      pp.get("thk",thk);
      buelerFlux->setViscosity(visc, thk, 2);
    }
    // TO BE IMPLEMENTED: GENERAL N LAYERS
    //else if ( nlayers>2 ) {
    //  Vector<Real> visc(nlayers); 
    //  pp.getarr("visc",visc,0,nlayers);
    //  Vector<Real> thk(nlayers-1);
    //  pp.getarr("thk",thk,0,nlayers-1);
    //}
    else { 
      MayDay::Error("Bueler flux nlayers not understood.");
    }
    
    Real flex;
    pp.get("flex", flex);
    buelerFlux->setFlexural(flex);

    bool includeElas = false;
    Real lame1=34.2666e9, lame2=26.6e9;
    pp.query("includeElas", includeElas);
    pp.query("lame1", lame1);
    pp.query("lame2", lame2);
    buelerFlux->setElastic(includeElas, lame1, lame2);

    bool oceanLoad = false;
    pp.query("oceanLoad", oceanLoad);
    buelerFlux->setOceanLoad(oceanLoad);

    bool ELRA = false;
    Real ELRAtau = 1e3; // in yr
    pp.query("ELRA", ELRA);
    pp.query("ELRAtau", ELRAtau);
    buelerFlux->setELRA(ELRA, ELRAtau);

    Real dt=0;
    pp.query("dt", dt);
    buelerFlux->setTimestep(dt);


    buelerFlux->precomputeGIAstep();

    // Check for initial disequilibrium uplift.
    // Assumes isostatic equilibrium is basal topography without ice
    // thickness. TO GENERALIZE AT A LATER DATE.
    bool init = false;
    pp.query("init", init);
    buelerFlux->setInitIceRef0( init );
    if (init) { 
      std::string initPrefix(a_prefix);
      initPrefix += ".init";
	  ParmParse pi(initPrefix.c_str());
      std::string file, name;
      pi.get("file", file);
      pi.get("name", name);
      // Read uplift
      RefCountedPtr<LevelData<FArrayBox>> upl0(new LevelData<FArrayBox>());
	  Vector<RefCountedPtr<LevelData<FArrayBox> > > data(1,upl0);
      Real dx;
      Vector<std::string> namevec(1,name);
	  readLevelData(data,dx,file,namevec,1);
      // Set uplift
      buelerFlux->setInitialVelocity(*upl0);
    }

    ptr = static_cast<SurfaceFlux*>(buelerFlux->new_surfaceFlux());

  }
#else
#warning ('buelerGIA flux not available')
  
#endif // BUELERGIA
#ifdef HAVE_PYTHON
  else if (type == "pythonFlux") {
    
    std::string module;
    pp.get("module",module);
    std::string function;
    pp.get("function",function);

    int nkwargs = 0;
    pp.query("n_kwargs",nkwargs);
    Vector<std::string> kwargName(nkwargs);
    if (nkwargs > 0)
      {
	pp.queryarr("kwargs",kwargName,0,nkwargs);
      }
    std::map<std::string, Real> kwarg;
    for (int i = 0; i < nkwargs; i++)
      {
	kwarg[kwargName[i]] = 0.0;
      }
    PythonInterface::PythonSurfaceFlux pythonFlux(module, function, kwarg);
    ptr = static_cast<SurfaceFlux*>(pythonFlux.new_surfaceFlux());

  }
#else
  #warning ('pythonFlux not available')
#endif
  else if (type == "")
    {
      ptr = NULL; // return a NULL and leave it up to the caller to care
    }
  else
    {
      // a type was specified but it made no sense...
      pout() << "unknown flux type " << type << std::endl;
      MayDay::Error("unknown flux type");
    }
  return ptr;
  
}


// default implementation -- in most cases this will return a Vector 
// of length 1, with a_root as the name
void SurfaceFlux::plot_names(const string& a_root, 
                             Vector<string>& a_plot_names) const
{
  int num_comps = num_plot_components();
  a_plot_names.resize(num_comps, a_root);
}

// default implementation simply calls evaluate function
void
SurfaceFlux::plot_data(LevelData<FArrayBox>& a_data,
                       const AmrIceBase& a_amrIce, 
                       int a_level, Real a_dt)
{
  evaluate(a_data, a_amrIce, a_level,a_dt);
}



#ifdef HAVE_PYTHON
#include "signal.h"


#endif
#include "NamespaceFooter.H"
