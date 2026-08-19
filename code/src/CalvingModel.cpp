#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "CalvingModel.H"
#include "MaskedCalvingModel.H"
#include "CrevasseCalvingModel.H"
#include "LevelMappedDerivatives.H"
#include "IceConstants.H"
#include "AmrIce.H"
#include "ParmParse.H"
#include "CalvingF_F.H"
#include "NamespaceHeader.H"

/// a default implementation
/**
   most models provide a criterion, rather than a rate.
 */
void
CalvingModel::getCalvingRate(LevelData<FArrayBox>& a_calvingRate, const AmrIce& a_amrIce,int a_level)
{
  DataIterator dit = a_calvingRate.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      a_calvingRate[dit].setVal(0.0);
    }
}

/// a default implementation
/**
   most models provide a criterion, rather than a rate.
*/ 
bool
CalvingModel::getCalvingVel(LevelData<FArrayBox>& a_centreCalvingVel,
			      const LevelData<FArrayBox>& a_centreIceVel,
			      const DisjointBoxLayout& a_grids,
			      const AmrIce& a_amrIce,int a_level)
{
  DataIterator dit = a_centreCalvingVel.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      a_centreCalvingVel[dit].setVal(0.0);
    }
  return true;
}

void
CalvingModel::getWaterDepth(LevelData<FArrayBox>& a_waterDepth, const AmrIce& a_amrIce,int a_level)
{
  DataIterator dit = a_waterDepth.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      a_waterDepth[dit].setVal(0.0);
    }
}

void
VariableRateCalvingModel::getCalvingRate(LevelData<FArrayBox>& a_calvingRate, const AmrIce& a_amrIce,int a_level)
{
  m_calvingRate->evaluate(a_calvingRate, a_amrIce, a_level, 0.0);
}


void 
DeglaciationCalvingModelA::applyCriterion
(LevelData<FArrayBox>& a_thickness, 
 LevelData<FArrayBox>& a_calvedIce,
 LevelData<FArrayBox>& a_addedIce,
 LevelData<FArrayBox>& a_removedIce, 
 LevelData<FArrayBox>& a_iceFrac, 
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{
  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
    {
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      FArrayBox& thck = a_thickness[dit];
      FArrayBox& calved = a_calvedIce[dit];
      FArrayBox& added = a_addedIce[dit];
      FArrayBox& removed = a_removedIce[dit];
      Box b = thck.box();

      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  Real prevThck = thck(iv);
	  if (mask(iv) == OPENSEAMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
	  else if (mask(iv) == OPENLANDMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
	  else
	    {
	      thck(iv) = std::max(thck(iv),m_minThickness);
	    }
	      
	  // Record gain/loss of ice
	  if (calved.box().contains(iv))
	    {
	      updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
	    }
	}
    }
}

void DomainEdgeCalvingModel::applyCriterion
(LevelData<FArrayBox>& a_thickness,
 LevelData<FArrayBox>& a_calvedIce,
 LevelData<FArrayBox>& a_addedIce,
 LevelData<FArrayBox>& a_removedIce,  
 LevelData<FArrayBox>& a_iceFrac, 
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{

  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  const DisjointBoxLayout& grids = levelCoords.grids();
  const ProblemDomain domain = grids.physDomain();
  const LevelData<BaseFab<int> >& levelMask = levelCoords.getFloatingMask();
  const IntVect ghost = a_thickness.ghostVect();
  //const LevelData<FArrayBox>& vt  = *a_amrIce.viscousTensor(a_level);
  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      //const Box& gridBox = grids[dit];
      for (int dir=0; dir<SpaceDim; dir++)
	{
	  if (!domain.isPeriodic(dir))
	    {

	      if (m_frontLo[dir] > 0)
		{
		  Box loBox = adjCellLo(domain,dir,ghost[dir]);
                  // (DFM 5-25-15) grow in transverse direction
                  // to ensure that we don't wind up with corner
                  // cells with ice in them
                  IntVect transverseVect = ghost;
                  transverseVect[dir] = 0;
                  loBox.grow(transverseVect);
		  loBox &= a_thickness[dit].box();
		  for (BoxIterator bit(loBox); bit.ok(); ++bit)
		    {
		      const IntVect& iv = bit();
		      const IntVect ip = iv + BASISV(dir);
		      //if (levelMask[dit](ip) != GROUNDEDMASKVAL)
		      Real prevThck = a_thickness[dit](iv);
		      a_thickness[dit](iv) = 0.0;
		      if (a_iceFrac[dit].box().contains(iv))
			a_iceFrac[dit](iv) = 0.0;
		      // Record gain/loss of ice
		      if (a_calvedIce[dit].box().contains(iv))
			{
			  updateCalvedIce(a_thickness[dit](iv),prevThck,levelMask[dit](iv),
					  a_addedIce[dit](iv),a_calvedIce[dit](iv),a_removedIce[dit](iv));
			}

		    }
		}
	      
	      if (m_frontHi[dir] > 0)
		{
		  Box hiBox = adjCellHi(domain,dir,ghost[dir]);
                  // (DFM 5-25-15) grow in transverse direction
                  // to ensure that we don't wind up with corner
                  // cells with ice in them
                  IntVect transverseVect = ghost;
                  transverseVect[dir] = 0;
                  hiBox.grow(transverseVect);
		  hiBox &= a_thickness[dit].box();
		  for (BoxIterator bit(hiBox); bit.ok(); ++bit)
		    {
		      const IntVect& iv = bit();
		      const IntVect ip = iv - BASISV(dir);
		      //if (levelMask[dit](ip) != GROUNDEDMASKVAL)
		      Real prevThck = a_thickness[dit](iv);
		      a_thickness[dit](iv) = 0.0;
		      if (a_iceFrac[dit].box().contains(iv))
			a_iceFrac[dit](iv) = 0.0;
		      // Record gain/loss of ice
		      if (a_calvedIce[dit].box().contains(iv))
			{
			  updateCalvedIce(a_thickness[dit](iv),prevThck,levelMask[dit](iv),
				      a_addedIce[dit](iv),a_calvedIce[dit](iv),a_removedIce[dit](iv));
			}

		    }
		} 
	    } // end if (!domain.isPeriodic(dir))
	} // end loop over dirs
      
      const BaseFab<int>& mask = levelMask[dit];
      FArrayBox& thck = a_thickness[dit];
      FArrayBox& calved = a_calvedIce[dit];
      FArrayBox& added = a_addedIce[dit];
      FArrayBox& removed = a_removedIce[dit];
      const Box& b = grids[dit];
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  Real prevThck = thck(iv);
	  if (m_preserveSea && mask(iv) == OPENSEAMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
	  else if (m_preserveLand && mask(iv) == OPENLANDMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
	  thck(iv) = std::max(thck(iv),0.0);

	  // Record gain/loss of ice
	  if (calved.box().contains(iv))
	    {
	      updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
	    }

	}

    } // end loop over boxes

}

void ProximityCalvingModel::applyCriterion
(LevelData<FArrayBox>& a_thickness,
 LevelData<FArrayBox>& a_calvedIce,
 LevelData<FArrayBox>& a_addedIce,
 LevelData<FArrayBox>& a_removedIce,  
 LevelData<FArrayBox>& a_iceFrac, 
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{

  Real time = a_amrIce.time();
  bool calvingActive = (time >= m_startTime && time < m_endTime);
  calvingActive = false;
  pout() << " time = " << time 
	 << " m_startTime = " <<  m_startTime
	 << " m_endTime = " <<  m_endTime
	 << "calvingActive = " << calvingActive
	 << std::endl;
  if (true || calvingActive)
    {
      const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
      const LevelData<FArrayBox>& proximity = *a_amrIce.groundingLineProximity(a_level);
      const LevelData<FArrayBox>& velocity = *a_amrIce.velocity(a_level);
      for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
	{
	  const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
	  FArrayBox& thck = a_thickness[dit];
	  FArrayBox& calved = a_calvedIce[dit];
	  FArrayBox& added = a_addedIce[dit];
	  FArrayBox& removed = a_removedIce[dit];
	  const FArrayBox& prox = proximity[dit];
	  const FArrayBox& vel = velocity[dit];
	  Box b = thck.box();b &= prox.box();
	  for (BoxIterator bit(b); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      Real prevThck = thck(iv);
	      Real vmod = std::sqrt(vel(iv,0)*vel(iv,0) + vel(iv,1)*vel(iv,1));
	      if (prox(iv) < m_proximity && calvingActive && vmod > m_velocity)
		{
		  //thck(iv) *= 0.5; thck(iv) = max(thck(iv),10.0);
		  thck(iv) = 0.0;
		}
	      if (mask(iv) == OPENSEAMASKVAL)
		{
		   thck(iv) = 0.0;
		}
	      if (mask(iv) == FLOATINGMASKVAL)
		{
		  thck(iv) = max(thck(iv),1.0);
		}

	      // Record gain/loss of ice
	      if (calved.box().contains(iv))
		{
		  updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
		}

	    }
	}
    }
}



CalvingModel* CalvingModel::parseCalvingModel(const char* a_prefix)
{

  CalvingModel* ptr = NULL;
  std::string type = "";
  ParmParse pp(a_prefix);
  pp.query("type",type);
  
  if (type == "NoCalvingModel")
    {
      ptr = new NoCalvingModel;
    }
  else if (type == "DomainEdgeCalvingModel")
    {
      Vector<int> frontLo(2,false); 
      pp.getarr("front_lo",frontLo,0,frontLo.size());
      Vector<int> frontHi(2,false);
      pp.getarr("front_hi",frontHi,0,frontHi.size());
      bool preserveSea = false;
      pp.query("preserveSea",preserveSea);
      bool preserveLand = false;
      pp.query("preserveLand",preserveLand);
      ptr = new DomainEdgeCalvingModel
	(frontLo, frontHi,preserveSea,preserveLand);
    }
  else if (type == "FixedFrontCalvingModel")
    {
      Real minThickness = 0.0;
      pp.get("min_thickness", minThickness );
      ptr = new DeglaciationCalvingModelA
	(0.0,  1.0e+10, minThickness, -1.2345678e+300, 1.2345678e+300);
    }
  else if (type == "DeglaciationCalvingModelA")
    {  
      Real minThickness = 0.0;
      pp.get("min_thickness", minThickness );
      Real calvingThickness = 0.0;
      pp.get("calving_thickness", calvingThickness );
      Real calvingDepth = 0.0;
      pp.get("calving_depth", calvingDepth );
      Real startTime = -1.2345678e+300;
      pp.query("start_time",  startTime);
      Real endTime = 1.2345678e+300;
      pp.query("end_time",  endTime);
      ptr = new DeglaciationCalvingModelA
	(calvingThickness,  calvingDepth, minThickness, startTime, endTime); 
    }
  else if (type == "DeglaciationCalvingModelB")
    {  
      Real minThickness = 0.0;
      pp.get("min_thickness", minThickness );
      Real calvingThickness = 0.0;
      pp.get("calving_thickness", calvingThickness );
      Real calvingDepth = 0.0;
      pp.query("calving_depth", calvingDepth );
      Real startTime = -1.2345678e+300;
      pp.query("start_time",  startTime);
      Real endTime = 1.2345678e+300;
      pp.query("end_time",  endTime);
      ptr = new DeglaciationCalvingModelB
	(calvingThickness,  calvingDepth, minThickness, startTime, endTime); 
    }
  else if (type == "ProximityCalvingModel")
    {
      Real proximity = 0.0;
      pp.get("proximity", proximity );
      Real velocity = 0.0;
      pp.query("velocity", velocity );
      Real startTime = -1.2345678e+300;
      pp.get("startTime",  startTime);
      Real endTime = 1.2345678e+300;
      pp.get("endTime",  endTime);
      ptr = new ProximityCalvingModel(proximity,velocity, startTime, endTime);
    }
  else if (type == "FlotationCalvingModel")
    {
      Vector<int> frontLo(2,false); 
      pp.getarr("front_lo",frontLo,0,frontLo.size());
      Vector<int> frontHi(2,false);
      pp.getarr("front_hi",frontHi,0,frontHi.size());
      bool preserveSea = false;
      pp.query("preserveSea",preserveSea);
      bool preserveLand = false;
      pp.query("preserveLand",preserveLand);
      ptr = new FlotationCalvingModel
	(frontLo, frontHi,preserveSea,preserveLand);
    }
  else if (type == "BennCalvingModel")
    {
      ptr = new BennCalvingModel(pp);
    }
  else if (type == "VanDerVeenCalvingModel")
    {
      ptr = new VdVCalvingModel(pp);
    }
  else if (type == "ThicknessCalvingModel")
    {  
      Real minThickness = 0.0;
      pp.get("min_thickness", minThickness );
      Real calvingThickness = 0.0;
      pp.get("calving_thickness", calvingThickness );
      Real calvingDepth = 0.0;
      pp.query("calving_depth", calvingDepth );
      Real startTime = -1.2345678e+300;
      pp.query("start_time",  startTime);
      Real endTime = 1.2345678e+300;
      pp.query("end_time",  endTime);
      bool factorMuCoef = false;
      pp.query("factor_mu_coef",factorMuCoef); 
      ptr = new ThicknessCalvingModel
	(calvingThickness,  calvingDepth, minThickness, startTime, endTime, factorMuCoef); 
    }
  else if (type == "CliffCollapseCalvingModel")
    {  
      Real maxCliffHeight = 100.0;
      pp.get("max_cliff_height", maxCliffHeight);
      Real recessionRate = 0.0;
      pp.get("recession_rate", recessionRate );
      Real startTime = -1.2345678e+300;
      pp.query("start_time",  startTime);
      Real endTime = 1.2345678e+300;
      pp.query("end_time",  endTime);
      ptr = new CliffCollapseCalvingModel
	(maxCliffHeight, recessionRate, startTime, endTime); 
    }
  else if (type == "MaxiumumExtentCalvingModel")  
    {
      Real startTime = -1.2345678e+300;
      pp.query("start_time",  startTime);
      Real endTime = 1.2345678e+300;
      pp.query("end_time",  endTime);

      Vector<Real> vect(SpaceDim,0.0);

      pp.getarr("lowLoc",vect,0,SpaceDim);
      RealVect lowLoc(D_DECL(vect[0], vect[1],vect[2]));      

      pp.getarr("highLoc",vect,0,SpaceDim);
      RealVect highLoc(D_DECL(vect[0], vect[1],vect[2]));      

      MaximumExtentCalvingModel* Ptr = new MaximumExtentCalvingModel(highLoc,
                                                                     lowLoc,
                                                                     startTime,
                                                                     endTime);
      ptr = static_cast<CalvingModel*>(Ptr);

    }
  else if (type == "MaskedCalvingModel")
    {
      Real minThickness = 0.0;
      pp.get("min_thickness", minThickness );

      // masked calving model uses a surfaceFlux as a mask
      std::string mask_prefix(a_prefix);
      mask_prefix += ".mask";
      SurfaceFlux* mask_ptr = SurfaceFlux::parse(mask_prefix.c_str());

      MaskedCalvingModel* Ptr = new MaskedCalvingModel(mask_ptr, minThickness);

      ptr = static_cast<CalvingModel*>(Ptr);

      // MaskedCalvingModel makes a copy of the mask, so clean up here
      if (mask_ptr != NULL)
        {
          delete mask_ptr;
        }
    }
  else if (type == "CompositeCalvingModel")
    {
      int nElements;
      pp.get("nElements",nElements);
     
      std::string elementPrefix(a_prefix);
      elementPrefix += ".element";

      Vector<CalvingModel*> elements(nElements);
      for (int i = 0; i < nElements; i++)
        {
          std::string prefix(elementPrefix);
          char s[32];
          sprintf(s,"%i",i);
          prefix += s;
          ParmParse pe(prefix.c_str());
          elements[i] = parseCalvingModel(prefix.c_str());
          CH_assert(elements[i] != NULL);
        }
      CompositeCalvingModel* compositePtr = new CompositeCalvingModel(elements);
      ptr = static_cast<CalvingModel*>(compositePtr);
    }
  
  else if (type == "VariableRateCalvingModel")
    {
      ptr = new VariableRateCalvingModel(pp);
    }

   else if (type == "RateAuBuhatCalvingModel")
    {
      ptr = new RateAuBuhatCalvingModel(pp);
    }
   else if (type == "RateProportionalToSpeedCalvingModel")
    {
  
      ptr = new RateAuBuhatCalvingModel(pp);
    }   
   else if (type == "VonMisesCalvingModel")
    {
      ptr = new VonMisesCalvingModel(pp);
    } 

  return ptr;
}


void 
DeglaciationCalvingModelB::applyCriterion
(LevelData<FArrayBox>& a_thickness,
 LevelData<FArrayBox>& a_calvedIce,
 LevelData<FArrayBox>& a_addedIce,
 LevelData<FArrayBox>& a_removedIce,  
 LevelData<FArrayBox>& a_iceFrac, 
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{
  
  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  
  for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
    {
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      FArrayBox& thck = a_thickness[dit];
      FArrayBox& calved = a_calvedIce[dit];
      FArrayBox& added = a_addedIce[dit];
      FArrayBox& removed = a_removedIce[dit];
      Box b = thck.box();
      
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  Real prevThck = thck(iv);
	  if (mask(iv) == OPENSEAMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
	  else if (mask(iv) == OPENLANDMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
          else if ((mask(iv) == FLOATINGMASKVAL) && (thck(iv) < m_calvingThickness))
            {
	      thck(iv) = m_minThickness;              
            }
	  else
	    {
	      thck(iv) = std::max(thck(iv),m_minThickness);
	    }

	  // Record gain/loss of ice
	  if (calved.box().contains(iv))
	    {
	      updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
	    }
	}
    }
}


void 
ThicknessCalvingModel::applyCriterion
(LevelData<FArrayBox>& a_thickness,
 LevelData<FArrayBox>& a_calvedIce,
 LevelData<FArrayBox>& a_addedIce,
 LevelData<FArrayBox>& a_removedIce,  
 LevelData<FArrayBox>& a_iceFrac, 
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{
  
  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
    {
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      FArrayBox& iceFrac = a_iceFrac[dit];
      FArrayBox& thck = a_thickness[dit];
      FArrayBox& calved = a_calvedIce[dit];
      FArrayBox& added = a_addedIce[dit];
      FArrayBox& removed = a_removedIce[dit];
      FArrayBox effectiveThickness(thck.box(), 1);
      effectiveThickness.copy(thck);

      if (m_factorMuCoef)
	{
	  effectiveThickness *= a_amrIce.muCoef(a_level)[dit];
	}

      
      Box b = thck.box();
      b &= iceFrac.box();
      
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();          
          // if iceFrac > 0, then rescale effectiveThickness
          // by dividing by iceFrac value, which gives "actual" thickness
          // in the partial cell. Probably eventually want to move this to 
          // fortran
	  Real prevThck = thck(iv);
          if (iceFrac(iv,0) > 0.0)
            {
              effectiveThickness(iv,0) /= iceFrac(iv,0);
            }
            
          if (mask(iv) == OPENLANDMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
          // allow ice to spread into open sea regions too, if appropriate
          else if (((mask(iv) == FLOATINGMASKVAL) || (mask(iv) == OPENSEAMASKVAL))
                   && (effectiveThickness(iv) < m_calvingThickness))
            {
              // note that we're setting thck here, not effectiveThickness, 
              // which is a temporary
              // also set the iceFrac to zero in these cells
	      thck(iv) = m_minThickness; 
              iceFrac(iv,0) = 0.0;
            }
	  else
	    {
	      thck(iv) = std::max(thck(iv),m_minThickness);
	    }

	  // Record gain/loss of ice
	  if (calved.box().contains(iv))
	    {
	      updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
	    }

	}
    }
}



  
//alter the thickness field at the end of a time step
void
MaximumExtentCalvingModel::applyCriterion(LevelData<FArrayBox>& a_thickness,
					  LevelData<FArrayBox>& a_calvedIce,
					  LevelData<FArrayBox>& a_addedIce,
					  LevelData<FArrayBox>& a_removedIce,  
					  LevelData<FArrayBox>& a_iceFrac, 
					  const AmrIce& a_amrIce,
					  int a_level,
					  Stage a_stage)
{

  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  const Real dx = a_amrIce.amrDx()[a_level];
  
  for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
    {
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      FArrayBox& thck = a_thickness[dit];
      FArrayBox& calved = a_calvedIce[dit];
      FArrayBox& added = a_addedIce[dit];
      FArrayBox& removed = a_removedIce[dit];
      Box b = thck.box();
      
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  Real prevThck = thck(iv);
          // compute location of cell center
          RealVect loc(iv);          
          loc += 0.5*RealVect::Unit;
          loc *= dx;
          
          // check high and low extents
          if ((mask(iv) == FLOATINGMASKVAL) || (mask(iv) == OPENSEAMASKVAL))
            {
              if (loc[0] <= m_lowLoc[0])
                {
                  thck(iv) = 0.0;
                }
              else if (loc[1] <= m_lowLoc[1])
                {
                  thck(iv) = 0.0;
                }
              else if (loc[0] > m_highLoc[0]) 
                {
                  thck(iv) = 0.0;
                }
              else if (loc[1] > m_highLoc[1])
                {
                  thck(iv) = 0.0;
                }
            } // end if floating or opensea

	  // Record gain/loss of ice
	  if (calved.box().contains(iv))
	    {
	      updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
	    }

        } // end loop over cells
  
    }

}




//alter the thickness field at the end of a time step
void 
CompositeCalvingModel::applyCriterion(LevelData<FArrayBox>& a_thickness, 
				      LevelData<FArrayBox>& a_calvedIce,
				      LevelData<FArrayBox>& a_addedIce,
				      LevelData<FArrayBox>& a_removedIce, 
				      LevelData<FArrayBox>& a_iceFrac, 
				      const AmrIce& a_amrIce,
				      int a_level,
				      Stage a_stage)
{
  for (int n=0; n<m_vectModels.size(); n++)
    {
      m_vectModels[n]->applyCriterion( a_thickness, a_calvedIce, a_addedIce, a_removedIce, a_iceFrac,a_amrIce, a_level, a_stage);
    }
}

//calving rate (scalar)
void 
CompositeCalvingModel::getCalvingRate(LevelData<FArrayBox>& a_rate, 
				      const AmrIce& a_amrIce,int a_level)
{

  m_vectModels[0]->getCalvingRate(a_rate, a_amrIce, a_level);
  LevelData<FArrayBox> tmp(a_rate.disjointBoxLayout(),1,a_rate.ghostVect());
  for (int n = 1; n < m_vectModels.size(); n++)
    {
      m_vectModels[n]->getCalvingRate(tmp, a_amrIce, a_level);
      for (DataIterator dit(a_rate.disjointBoxLayout()); dit.ok(); ++dit)
	{
	  a_rate[dit] += tmp[dit];
	}
    }
}

bool 
CompositeCalvingModel::getCalvingVel
(LevelData<FArrayBox>& a_centreCalvingVel,
 const LevelData<FArrayBox>& a_centreIceVel,
 const DisjointBoxLayout& a_grids,
 const AmrIce& a_amrIce,int a_level)
{
  bool s = m_vectModels[0]->getCalvingVel(a_centreCalvingVel, a_centreIceVel,
					  a_grids, a_amrIce, a_level);
  int n = 1;
  LevelData<FArrayBox> tmp(a_grids,SpaceDim,a_centreCalvingVel.ghostVect());
  while ((s) && (n <  m_vectModels.size()))
  {
    s = m_vectModels[n]->getCalvingVel(tmp, a_centreIceVel,
				       a_grids, a_amrIce, a_level);
    if (s)
      {
	for (DataIterator dit(a_grids); dit.ok(); ++dit)
	  {
	    a_centreCalvingVel[dit] += tmp[dit];
	  }
      }
    n++;
  }
  return s;
}
 
CompositeCalvingModel::~CompositeCalvingModel()
{
  for (int n=0; n<m_vectModels.size(); n++)
    {
      delete m_vectModels[n];
      m_vectModels[n] = NULL;
    }
}

void FlotationCalvingModel::applyCriterion
(LevelData<FArrayBox>& a_thickness,
 LevelData<FArrayBox>& a_calvedIce,
 LevelData<FArrayBox>& a_addedIce,
 LevelData<FArrayBox>& a_removedIce,  
 LevelData<FArrayBox>& a_iceFrac, 
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{

  m_domainEdgeCalvingModel.applyCriterion( a_thickness, a_calvedIce, a_addedIce, a_removedIce, a_iceFrac,a_amrIce, a_level, a_stage);
  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
    {
      FArrayBox& thck = a_thickness[dit];
      FArrayBox& calved = a_calvedIce[dit];
      FArrayBox& added = a_addedIce[dit];
      FArrayBox& removed = a_removedIce[dit];
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      const Box& b = levelCoords.grids()[dit];
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  Real prevThck = thck(iv);
	  if (mask(iv) == FLOATINGMASKVAL)
	    {
	      thck(iv) = 0.0; 
	    }

	  // Record gain/loss of ice
	  if (calved.box().contains(iv))
	    {
	      updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
	    }

	}
    }
}

void
CalvingModel::updateCalvedIce(const Real& a_thck, const Real a_prevThck, const int a_mask, Real& a_added, Real& a_calved, Real& a_removed)
{

  if (a_thck > a_prevThck)
    {
      a_added += (a_prevThck-a_thck);
    }
  else 
    {
      if ((a_mask == OPENSEAMASKVAL) || (a_mask == FLOATINGMASKVAL))
	{
	  a_calved += (a_prevThck-a_thck);
	}
      else
	{
	  a_calved += (a_prevThck-a_thck);
	}
    } 

}


VariableRateCalvingModel::VariableRateCalvingModel(ParmParse& a_pp)
{
      Real startTime = -1.2345678e+300;
      a_pp.query("start_time",  startTime);
      Real endTime = 1.2345678e+300;
      a_pp.query("end_time",  endTime);
 
      Vector<int> frontLo(2,false); 
      a_pp.getarr("front_lo",frontLo,0,frontLo.size());
      Vector<int> frontHi(2,false);
      a_pp.getarr("front_hi",frontHi,0,frontHi.size());
      bool preserveSea = false;
      a_pp.query("preserveSea",preserveSea);
      bool preserveLand = false;
      a_pp.query("preserveLand",preserveLand);

      m_domainEdgeCalvingModel = new DomainEdgeCalvingModel(frontLo,frontHi,preserveSea,preserveLand);

      std::string prefix (a_pp.prefix());
      m_calvingRate = SurfaceFlux::parse( (prefix + ".CalvingRate").c_str());

}

void VariableRateCalvingModel::applyCriterion
(LevelData<FArrayBox>& a_thickness,
 LevelData<FArrayBox>& a_calvedIce,
 LevelData<FArrayBox>& a_addedIce,
 LevelData<FArrayBox>& a_removedIce,  
 LevelData<FArrayBox>& a_iceFrac, 
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{

  (*m_domainEdgeCalvingModel).applyCriterion( a_thickness, a_calvedIce, a_addedIce, a_removedIce, a_iceFrac,a_amrIce, a_level, a_stage);

  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
    {
      FArrayBox& thck = a_thickness[dit];
      FArrayBox& calved = a_calvedIce[dit];
      FArrayBox& added = a_addedIce[dit];
      FArrayBox& removed = a_removedIce[dit];
      FArrayBox& frac = a_iceFrac[dit];
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      Real frac_eps = TINY_FRAC;
      const Box& b = levelCoords.grids()[dit];
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  Real prevThck = thck(iv);

	  if (frac(iv) < frac_eps * frac_eps)
	    {
	      frac(iv) = 0.0;
	    }
	  
	  if (frac(iv) < frac_eps)
	    {
	      thck(iv)=0.0;
	    }

	  // Record gain/loss of ice
	  if (calved.box().contains(iv))
	    {
	      updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
	    }

	}
    }
}



CalvingModel* VariableRateCalvingModel::new_CalvingModel()
  {
    VariableRateCalvingModel* ptr = new VariableRateCalvingModel(*this);
    ptr->m_startTime = m_startTime;
    ptr->m_endTime = m_endTime;
    ptr->m_calvingRate = m_calvingRate->new_surfaceFlux();
    ptr->m_domainEdgeCalvingModel = new DomainEdgeCalvingModel(*m_domainEdgeCalvingModel);
    return ptr; 
  }

VariableRateCalvingModel::~VariableRateCalvingModel()
{

  if (m_domainEdgeCalvingModel != NULL)
    {
      delete m_domainEdgeCalvingModel;
      m_domainEdgeCalvingModel = NULL;
    }

  if (m_calvingRate != NULL)
    {
      delete m_calvingRate;
      m_calvingRate = NULL;
    }

}


void 
CliffCollapseCalvingModel::applyCriterion(LevelData<FArrayBox>& a_thickness,
					  LevelData<FArrayBox>& a_calvedIce,
					  LevelData<FArrayBox>& a_addedIce,
					  LevelData<FArrayBox>& a_removedIce,  
					  LevelData<FArrayBox>& a_iceFrac, 
					  const AmrIce& a_amrIce,
					  int a_level,
					  Stage a_stage)
{

  // only do this at the end of a timestep
  // (since that's the only time a time-integrated recession rate makes any sense)
  if (a_stage == PostThicknessAdvection)
    {
      Real dt = a_amrIce.dt();
      Real dx = a_amrIce.amrDx()[a_level];
      
      const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
      const LevelData<FArrayBox>& surfaceHeight = levelCoords.getSurfaceHeight();
      
      for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
	{
	  const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
	  FArrayBox& iceFrac = a_iceFrac[dit];
	  FArrayBox& thck = a_thickness[dit];
	  FArrayBox& calved = a_calvedIce[dit];
	  FArrayBox& added = a_addedIce[dit];
	  FArrayBox& removed = a_removedIce[dit];
	  const FArrayBox& surface = surfaceHeight[dit];
	  FArrayBox effectiveSurface(surface.box(),1);
	  effectiveSurface.copy(surface);
	  FArrayBox effectiveThickness(thck.box(), 1);
	  effectiveThickness.copy(thck);
	  Box b = thck.box();
	  b &= iceFrac.box();
	  
	  // keep track of which cells we've already done in order to avoid double-counting
	  BaseFab<int> alreadyDone(b,1);
	  alreadyDone.setVal(0);
	  
	  Real phiNew;
	  
	  for (BoxIterator bit(b); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();          
	      // if iceFrac > 0, then rescale effectiveThickness
	      // by dividing by iceFrac value, which gives "actual" thickness
	      // in the partial cell. Probably eventually want to move this to 
	      // fortran
	      // also compute "effective surface", which is the upper surface height based
	      // on the effective thickness rather than the cell-averaged thickness
	      // Probably eventually want to move this to fortran
	      Real prevThck = thck(iv);
	      if (iceFrac(iv,0) > 0.0)
		{
		  effectiveThickness(iv,0) /= iceFrac(iv,0);
		  effectiveSurface(iv,0) += effectiveThickness(iv,0) - thck(iv,0);	      
		}
	      
	      // if ice is grounded, look at neighbors to see if there are any empty neighbors, then look at
	      // surface differences.
	      if (mask(iv) == GROUNDEDMASKVAL) 
		{
		  // loop over directions
		  for (int dir=0; dir<SpaceDim; dir++)
		    {
		      IntVect shiftVect = BASISV(dir);
		      IntVect ivp = iv + shiftVect;
		      IntVect ivm = iv - shiftVect;
		      
		      // look in both high and low directions at once
		      if (((mask(ivp,0) != GROUNDEDMASKVAL) && (mask(ivp,0) != OPENLANDMASKVAL) &&  ((effectiveSurface(iv,0) - effectiveSurface(ivp,0)) > m_maxCliffThickness)) ||
			  ((mask(ivm,0) != GROUNDEDMASKVAL) && (mask(ivm,0) != OPENLANDMASKVAL) && ((effectiveSurface(iv,0) - effectiveSurface(ivm,0)) > m_maxCliffThickness)))
			{
			  // we have a cliff!  only adjust this cell if we haven't already
			  if (alreadyDone(iv,0) == 0)
			    {
			      alreadyDone(iv,0) = 1;
			      phiNew = iceFrac(iv,0) - m_recessionRate*dt/dx;
			      // don't go below zero
			      phiNew = Max(phiNew, 0.0);
			      
			      // note that we're setting thck here, not effectiveThickness, 
			      // which is a temporary
			      // also modify the iceMask to zero in these cells
			      
			      thck(iv,0) = thck(iv,0)*phiNew/iceFrac(iv,0);
			      iceFrac(iv,0) = phiNew;
			    } // end we haven't already done this one
			} // end if we have a cliff
		    } // end loop over directions
		} // end if this cell is grounded (no floating cliffs)
	      
	      
	      // Record gain/loss of ice
	      if (calved.box().contains(iv))
		{
		  updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
		}
	      
	    } // end loop over cells in this box
	} // end loop over grids on this level	  
    } // end if we're at the post-advection stage
}


RateAuBuhatCalvingModel::RateAuBuhatCalvingModel(ParmParse& a_pp)
{
      Real startTime = -1.2345678e+300;
      a_pp.query("start_time",  startTime);
      Real endTime = 1.2345678e+300;
      a_pp.query("end_time",  endTime);
 
      Vector<int> frontLo(2,false); 
      a_pp.getarr("front_lo",frontLo,0,frontLo.size());
      Vector<int> frontHi(2,false);
      a_pp.getarr("front_hi",frontHi,0,frontHi.size());
      bool preserveSea = false;
      a_pp.query("preserveSea",preserveSea);
      bool preserveLand = false;
      a_pp.query("preserveLand",preserveLand);

      m_domainEdgeCalvingModel = new DomainEdgeCalvingModel(frontLo,frontHi,preserveSea,preserveLand);

      std::string prefix (a_pp.prefix());
      m_proportion = SurfaceFlux::parse( (prefix + ".proportion").c_str());
      if (!m_proportion) m_proportion = new zeroFlux(); 
      m_independent = SurfaceFlux::parse( (prefix + ".independent").c_str());
      if (!m_independent) m_independent  = new zeroFlux();
      m_independent_normal = false;
      a_pp.query("independent_normal",m_independent_normal); 
      m_vector = true; // this is essentially required
      a_pp.query("vector", m_vector);

      
}

void RateAuBuhatCalvingModel::applyCriterion
(LevelData<FArrayBox>& a_thickness,
 LevelData<FArrayBox>& a_calvedIce,
 LevelData<FArrayBox>& a_addedIce,
 LevelData<FArrayBox>& a_removedIce,  
 LevelData<FArrayBox>& a_iceFrac, 
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{
  // No explicit criterion in this case, but m_domainEdgeCalvingModel applies.
  (*m_domainEdgeCalvingModel).applyCriterion( a_thickness, a_calvedIce, a_addedIce, a_removedIce, a_iceFrac,a_amrIce, a_level, a_stage);
}



CalvingModel* RateAuBuhatCalvingModel::new_CalvingModel()
  {
    RateAuBuhatCalvingModel* ptr = new RateAuBuhatCalvingModel(*this);
    ptr->m_startTime = m_startTime;
    ptr->m_endTime = m_endTime;
    ptr->m_proportion = m_proportion->new_surfaceFlux();
    ptr->m_independent = m_independent->new_surfaceFlux();
    ptr->m_domainEdgeCalvingModel = new DomainEdgeCalvingModel(*m_domainEdgeCalvingModel);
    return ptr; 
  }

RateAuBuhatCalvingModel::~RateAuBuhatCalvingModel()
{

  if (m_domainEdgeCalvingModel != NULL)
    {
      delete m_domainEdgeCalvingModel;
      m_domainEdgeCalvingModel = NULL;
    }

  if (m_proportion != NULL)
    {
      delete m_proportion;
      m_proportion = NULL;
    }

    if (m_independent != NULL)
    {
      delete m_independent;
      m_independent = NULL;
    }

  
}

bool 
RateAuBuhatCalvingModel::getCalvingVel
(LevelData<FArrayBox>& a_centreCalvingVel,
 const LevelData<FArrayBox>& a_centreIceVel,
 const DisjointBoxLayout& a_grids,
 const AmrIce& a_amrIce,int a_level)
{
  if (!m_vector) return false;
  
  // cell-centered proportion
  LevelData<FArrayBox> prop(a_grids, 1, 2*IntVect::Unit);
  m_proportion->evaluate(prop, a_amrIce, a_level, a_amrIce.dt());
  prop.exchange();
  
  // -velocity * proportion
  for (DataIterator dit(a_grids); dit.ok(); ++dit)
    {
      prop[dit] *= -1;
      a_centreCalvingVel[dit].copy(a_centreIceVel[dit]);
      for (int dir = 0; dir < SpaceDim; ++dir)
	{
	  a_centreCalvingVel[dit].mult(prop[dit], 0, dir, 1);
	}
    }

  if (m_independent)
    {
      // cell-centered independent part 
      LevelData<FArrayBox> ccrate(a_grids, 1, 1*IntVect::Unit);
      m_independent->evaluate(ccrate, a_amrIce, a_level, a_amrIce.dt());
      ccrate.exchange();
     
      for (DataIterator dit(a_grids); dit.ok(); ++dit)
        {
	  const FArrayBox& u = a_centreIceVel[dit];
	  const FArrayBox& frac = (*a_amrIce.iceFraction(a_level))[dit];
	  FArrayBox& u_c = a_centreCalvingVel[dit];
	  
	  Box gbox = a_grids[dit];
	  gbox.grow(1);
	  for (BoxIterator bit(gbox); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      if (m_independent_normal)
	      {
		Real dfdx = frac(iv + BASISV(0)) - frac(iv - BASISV(0));
		Real dfdy = frac(iv + BASISV(1)) - frac(iv - BASISV(1));
		Real df = 1.0e-10 + std::sqrt(dfdx*dfdx + dfdy*dfdy);
		u_c(iv,0) += ccrate[dit](iv)*dfdx/df;
		u_c(iv,1) += ccrate[dit](iv)*dfdy/df;
	      }
	      else
	      {	      
	      	Real umod = 1.0e-10 + std::sqrt(u(iv,0)*u(iv,0) + u(iv,1)*u(iv,1));
	      	u_c(iv,0) -=  ccrate[dit](iv)*u(iv,0) / umod;
	      	u_c(iv,1) -=  ccrate[dit](iv)*u(iv,1) / umod;
	      } // normal vs anti-parallel
	    } // bit
	} // dit
    } // m_independent
  return true;
  
}



void
RateAuBuhatCalvingModel::getCalvingRate
(LevelData<FArrayBox>& a_calvingRate, const AmrIce& a_amrIce,int a_level)
{
  // CH_assert(false); // we don't want to use this
  m_proportion->evaluate(a_calvingRate, a_amrIce, a_level, a_amrIce.dt());
  LevelData<FArrayBox> indep(a_calvingRate.disjointBoxLayout(), 1, 2*IntVect::Unit);
  if (m_independent) m_independent->evaluate(indep, a_amrIce, a_level, a_amrIce.dt());
  const LevelData<FArrayBox>& vel = *a_amrIce.velocity(a_level); // flux vel might be better  
  for (DataIterator dit(vel.dataIterator()); dit.ok(); ++dit)
    {
      Box b = a_calvingRate[dit].box();
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  Real usq = 0.0;
	  for (int dir = 0; dir < SpaceDim; dir++)
	    {
	      usq += vel[dit](iv,dir)*vel[dit](iv,dir);
	    }	
	  a_calvingRate[dit](iv) *= std::sqrt(usq);
	}
      if  (m_independent)
	{
	  a_calvingRate[dit] += indep[dit];
	}
    }

  int dbg = 0;dbg++; 
}

VonMisesCalvingModel::VonMisesCalvingModel(ParmParse& a_pp)
{
      Real startTime = -1.2345678e+300;
      a_pp.query("start_time",  startTime);
      Real endTime = 1.2345678e+300;
      a_pp.query("end_time",  endTime);
 
      Vector<int> frontLo(2,false); 
      a_pp.getarr("front_lo",frontLo,0,frontLo.size());
      Vector<int> frontHi(2,false);
      a_pp.getarr("front_hi",frontHi,0,frontHi.size());
      bool preserveSea = false;
      a_pp.query("preserveSea",preserveSea);
      bool preserveLand = false;
      a_pp.query("preserveLand",preserveLand);

      m_domainEdgeCalvingModel = new DomainEdgeCalvingModel(frontLo,frontHi,preserveSea,preserveLand);

      std::string prefix (a_pp.prefix());
      m_scale = SurfaceFlux::parse( (prefix + ".scale").c_str());
      if (!m_scale) m_scale = new zeroFlux(); 
      m_independent = SurfaceFlux::parse( (prefix + ".independent").c_str());
      if (!m_independent) m_independent  = new zeroFlux(); 
      m_vector = true;
      a_pp.query("vector", m_vector);

      
}

CalvingModel* VonMisesCalvingModel::new_CalvingModel()
  {
    VonMisesCalvingModel* ptr = new VonMisesCalvingModel(*this);
    ptr->m_startTime = m_startTime;
    ptr->m_endTime = m_endTime;
    ptr->m_scale = m_scale->new_surfaceFlux();
    ptr->m_independent = m_independent->new_surfaceFlux();
    ptr->m_domainEdgeCalvingModel = new DomainEdgeCalvingModel(*m_domainEdgeCalvingModel);
    return ptr; 
  }

VonMisesCalvingModel::~VonMisesCalvingModel()
{

  if (m_domainEdgeCalvingModel != NULL)
    {
      delete m_domainEdgeCalvingModel;
      m_domainEdgeCalvingModel = NULL;
    }

  if (m_scale != NULL)
    {
      delete m_scale;
      m_scale = NULL;
    }

    if (m_independent != NULL)
    {
      delete m_independent;
      m_independent = NULL;
    }

  
}

	
void VonMisesCalvingModel::applyCriterion
(LevelData<FArrayBox>& a_thickness,
 LevelData<FArrayBox>& a_calvedIce,
 LevelData<FArrayBox>& a_addedIce,
 LevelData<FArrayBox>& a_removedIce,  
 LevelData<FArrayBox>& a_iceFrac, 
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{
  // No explicit criterion in this case, but m_domainEdgeCalvingModel applies.
  (*m_domainEdgeCalvingModel).applyCriterion( a_thickness, a_calvedIce, a_addedIce, a_removedIce, a_iceFrac,a_amrIce, a_level, a_stage);
}


bool 
VonMisesCalvingModel::getCalvingVel
(LevelData<FArrayBox>& a_centreCalvingVel,
 const LevelData<FArrayBox>& a_centreIceVel,
 const DisjointBoxLayout& a_grids,
 const AmrIce& a_amrIce,int a_level)
{
  if (!m_vector) return false;
  
  // cell-centered scale
  LevelData<FArrayBox> scale(a_grids, 1, 2*IntVect::Unit);
  m_scale->evaluate(scale, a_amrIce, a_level, a_amrIce.dt());
  scale.exchange();

  LevelData<FArrayBox> vonmises(a_grids, 1, 2*IntVect::Unit);
  const LevelData<FArrayBox>& viscousTensor = *a_amrIce.viscousTensor(a_level);
  const LevelSigmaCS& geometry = *a_amrIce.geometry(a_level);
  const LevelData<FArrayBox>& thickness = geometry.getH();
  
  // locate specific components of the viscous tensor the multicomponent arrays.
  // the derivComponent function is in LevelMappedDerivatives.
  int xxComp = derivComponent(0,0);
  int xyComp = derivComponent(1,0);
  int yxComp = derivComponent(0,1);
  int yyComp = derivComponent(1,1);

  Real eps = 1.0e-10;
  
  // Loop over individual boxes on a single processor
  DataIterator dit=a_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      // valid region box for this patch
      const Box& thisBox = a_grids[dit];
      //vonmises[dit].setVal(0.0); 
      // Compute the Von Mises Stress (cell-centered)
      FORT_VONMISES(CHF_FRA1(vonmises[dit],0),
		    CHF_CONST_FRA(viscousTensor[dit]),
		    CHF_CONST_FRA1(thickness[dit],0),
		    CHF_INT(xxComp),
		    CHF_INT(xyComp),
		    CHF_INT(yxComp),
		    CHF_INT(yyComp),
		    CHF_CONST_REAL(eps),
		    CHF_BOX(thisBox));
      
    } // End loop over boxes on a single processor
  // End Von Mises calculation

  // -velocity * sigma_vm * scale
  for (DataIterator dit(a_grids); dit.ok(); ++dit)
    {
      const FArrayBox& u = a_centreIceVel[dit];
      FArrayBox& v = a_centreCalvingVel[dit];
      v.copy(u);
      for (int dir = 0; dir < SpaceDim; ++dir)
	{
	  v.mult(vonmises[dit], 0, dir, 1);
	  v.mult(scale[dit], 0, dir, 1);
	  
	  // \todo fix this, just a test. avoid -v > 2u
	  for (BoxIterator bit(v.box()); bit.ok(); ++bit)
	  {
		const IntVect& iv = bit();
		if (Abs(v(iv,dir)) > 2.0*Abs(u(iv,dir)))
		{
			v(iv,dir) = 2.0*u(iv,dir);
		}
	  }	  
	}	
	v *= -1; // opposing direction
      
				
    }
  
  return true;
  
}


void
VonMisesCalvingModel::getCalvingRate
(LevelData<FArrayBox>& a_calvingRate, const AmrIce& a_amrIce,int a_level)
{
  const DisjointBoxLayout& a_grids = a_calvingRate.disjointBoxLayout();
  // CH_assert(false); // we don't want to use this
  m_scale->evaluate(a_calvingRate, a_amrIce, a_level, a_amrIce.dt());
  LevelData<FArrayBox> indep(a_calvingRate.disjointBoxLayout(), 1, 2*IntVect::Unit);
  if (m_independent) m_independent->evaluate(indep, a_amrIce, a_level, a_amrIce.dt());

  // Von Mises Block
  LevelData<FArrayBox> vonmises(a_calvingRate.disjointBoxLayout(), 1, 2*IntVect::Unit);
  // Access the (cell-centered) viscous tensor (vertically integrated stress)
  // SHOULD USE FACE VISCOUS TENSOR!
  const LevelData<FArrayBox>& a_viscousTensor = *a_amrIce.viscousTensor(a_level);
  //const LevelData<FArrayBox> &a_thickness = (*a_amrIce.geometry(a_level)).a_coordSys.getH();
  const LevelSigmaCS& a_geometry = *a_amrIce.geometry(a_level);
  const LevelData<FArrayBox>& a_thickness = a_geometry.getH();

  // locate specific components of the viscous tensor the multicomponent arrays.
  // the derivComponent function is in LevelMappedDerivatives.
  int xxComp = derivComponent(0,0);
  int xyComp = derivComponent(1,0);
  int yxComp = derivComponent(0,1);
  int yyComp = derivComponent(1,1);

  Real eps = 1.0e-10;

  // Loop over individual boxes on a single processor
  DataIterator dit=a_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      // valid region box for this patch
      const Box& thisBox = a_grids[dit];

  // Compute the Von Mises Stress (cell-centered)
  FORT_VONMISES(CHF_FRA1(vonmises[dit],0),
          CHF_CONST_FRA(a_viscousTensor[dit]),
          CHF_CONST_FRA1(a_thickness[dit],0),
          CHF_INT(xxComp),
          CHF_INT(xyComp),
          CHF_INT(yxComp),
          CHF_INT(yyComp),
          CHF_CONST_REAL(eps),
          CHF_BOX(thisBox));

    } // End loop over boxes on a single processor
  // End Von Mises block


  const LevelData<FArrayBox>& vel = *a_amrIce.velocity(a_level); // flux vel might be better  
  for (DataIterator dit(vel.dataIterator()); dit.ok(); ++dit)
    {
      Box b = a_calvingRate[dit].box();
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  Real usq = 0.0;
	  for (int dir = 0; dir < SpaceDim; dir++)
	    {
	      usq += vel[dit](iv,dir)*vel[dit](iv,dir);
	    }	
	  a_calvingRate[dit](iv) *= std::sqrt(usq)*vonmises[dit](iv);
	}
      if  (m_independent)
	{
	  a_calvingRate[dit] += indep[dit];
	}
    }

  int dbg = 0;dbg++; 

}
/**Von Mises Building Blocks**/

#include "NamespaceFooter.H"
