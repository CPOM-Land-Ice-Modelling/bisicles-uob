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


void DomainEdgeCalvingModel::evaluateCriterion
(LevelData<BaseFab<bool > >& a_critical,
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{

  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  const DisjointBoxLayout& grids = levelCoords.grids();
  const ProblemDomain domain = grids.physDomain();
  const LevelData<BaseFab<int> >& levelMask = levelCoords.getFloatingMask();
  const IntVect ghost = a_critical.ghostVect();
  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
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
		  loBox &= a_critical[dit].box();
		  for (BoxIterator bit(loBox); bit.ok(); ++bit)
		    {
		      const IntVect& iv = bit();
		      a_critical[dit](iv) = true;
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
		  hiBox &= a_critical[dit].box();
		  for (BoxIterator bit(hiBox); bit.ok(); ++bit)
		    {
		      const IntVect& iv = bit();
		      a_critical[dit](iv) = true;
		    }
		} 
	    } // end if (!domain.isPeriodic(dir))
	} // end loop over dirs
      
      const BaseFab<int>& mask = levelMask[dit];
      const Box& b = grids[dit];
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  if (m_preserveSea && mask(iv) == OPENSEAMASKVAL)
	    {
	      a_critical[dit](iv) = true;
	    }
	  else if (m_preserveLand && mask(iv) == OPENLANDMASKVAL)
	    {
	      a_critical[dit](iv) = true;
	    }
	}

    } // end loop over boxes

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
      ptr = new ThicknessCalvingModel
	(0.0,  1.0e+10, minThickness, -1.2345678e+300, 1.2345678e+300, false);
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
ThicknessCalvingModel::evaluateCriterion
(LevelData<BaseFab<bool > >& a_critical,
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{
  
  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  const LevelData<FArrayBox>& iceFracData = *a_amrIce.iceFrac(a_level);
  const LevelData<FArrayBox>& thckData = levelCoords.getH();
  for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
    {
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      const FArrayBox& iceFrac = iceFracData[dit];
      const FArrayBox& thck = thckData[dit];
      FArrayBox effectiveThickness(thck.box(), 1);
      effectiveThickness.copy(thck);
      if (m_factorMuCoef)
	{
	  effectiveThickness *= a_amrIce.muCoef(a_level)[dit];
	}
      Box b = thck.box();
      b &= iceFrac.box();
      b &= a_critical[dit].box();
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();          
	  if (iceFrac(iv,0) > 0.0)
	    {
	      effectiveThickness(iv,0) /= iceFrac(iv,0);
	    }
	  if (mask(iv) == OPENLANDMASKVAL)
	    {
	      a_critical[dit](iv) = true;
	    }
          else if (((mask(iv) == FLOATINGMASKVAL) || (mask(iv) == OPENSEAMASKVAL))
                   && (effectiveThickness(iv) < m_calvingThickness))
            {
	      a_critical[dit](iv) = true;
            }
	}
    }
}



  
//alter the thickness field at the end of a time step
void 
CompositeCalvingModel::evaluateCriterion(LevelData<BaseFab<bool > >& a_critical,
				      const AmrIce& a_amrIce,
				      int a_level,
				      Stage a_stage)
{
  for (int n=0; n<m_vectModels.size(); n++)
    {
      LevelData<BaseFab<bool> > critical(a_critical.disjointBoxLayout(),1,a_critical.ghostVect());

      m_vectModels[n]->evaluateCriterion( critical, a_amrIce, a_level, a_stage);
      for (DataIterator dit(critical.disjointBoxLayout()); dit.ok(); ++dit)
	{
	  BaseFab<bool>& crit  = critical[dit];
	  BaseFab<bool>& critAcc = a_critical[dit];
	  const Box& b = critAcc.box();
	  for (BoxIterator bit(b); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      critAcc(iv) = critAcc(iv) || crit(iv);
	    }
	}
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

void FlotationCalvingModel::evaluateCriterion
(LevelData<BaseFab<bool > >& a_critical,
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{

  m_domainEdgeCalvingModel.evaluateCriterion( a_critical, a_amrIce, a_level, a_stage);
  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
    {
      BaseFab<bool>& crit = a_critical[dit];
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      const Box& b = levelCoords.grids()[dit];
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  if (mask(iv) == FLOATINGMASKVAL)
	    {
	      crit(iv) = true; 
	    }
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

void VariableRateCalvingModel::evaluateCriterion
(LevelData<BaseFab<bool > >& a_critical,
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{

  (*m_domainEdgeCalvingModel).evaluateCriterion( a_critical, a_amrIce, a_level, a_stage);

  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  const LevelData<FArrayBox>& iceFracData = *a_amrIce.iceFrac(a_level);
  for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
    {
      const FArrayBox& frac = iceFracData[dit];
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      BaseFab<bool>& crit = a_critical[dit];
      Real frac_eps = TINY_FRAC;
      const Box& b = levelCoords.grids()[dit];
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  if (frac(iv) < frac_eps)
	    {
	      crit(iv) = true;
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

void RateAuBuhatCalvingModel::evaluateCriterion
(LevelData<BaseFab<bool > >& a_critical,
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{
  // No explicit criterion in this case, but m_domainEdgeCalvingModel applies.
  (*m_domainEdgeCalvingModel).evaluateCriterion( a_critical, a_amrIce, a_level, a_stage);
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

	
void VonMisesCalvingModel::evaluateCriterion
(LevelData<BaseFab<bool > >& a_critical,
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{
  // No explicit criterion in this case, but m_domainEdgeCalvingModel applies.
  (*m_domainEdgeCalvingModel).evaluateCriterion( a_critical, a_amrIce, a_level, a_stage);
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
