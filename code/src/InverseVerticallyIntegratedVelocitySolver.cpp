#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "InverseVerticallyIntegratedVelocitySolver.H"
#include "CGOptimize.H"
#include "ParmParse.H"
#include "ControlF_F.H"
#include "BisiclesF_F.H"
#include "ReflectGhostCells.H"
#include "QuadCFInterp.H"
#include "PiecewiseLinearFillPatch.H"
#include "CoarseAverage.H"
#include "CoarseAverageFace.H"
#include "IceUtility.H"
#include "IceConstants.H"
#include "computeNorm.H"
#include "computeSum.H"
#include "LevelDataBasalFriction.H"
#include "PatchGodunov.H"
#include "AdvectPhysics.H"
#include "DivergenceF_F.H"
#include "AmrIceBase.H"
#include "NamespaceHeader.H"

// for writing l-curve params to csv file
#include <iostream>
#include <fstream>
#include <vector>

#define CCOMP 0
#define MUCOMP 1
#define NXCOMP MUCOMP + 1



InverseVerticallyIntegratedVelocitySolver::Configuration::Configuration()
  :m_velObs_c(NULL),
   m_velObs_x(NULL),
   m_damRegCoef_c(NULL),
   m_X1TimeRegCoef_c(NULL),
#if CH_SPACEDIM > 1
   m_velObs_y(NULL), 
#endif
   m_gradientFactor(NULL)
{

}

InverseVerticallyIntegratedVelocitySolver::Configuration::~Configuration()
{

  if (m_velObs_x != NULL)
    {
      delete m_velObs_x; m_velObs_x = NULL;
    }
#if CH_SPACEDIM > 1
  if (m_velObs_y != NULL)
    {
      delete m_velObs_y; m_velObs_y = NULL;
    }
#endif
  if (m_velObs_c != NULL)
    {
      delete m_velObs_c; m_velObs_c = NULL;
    }
    
    if (m_damRegCoef_c != NULL)
    {
      delete m_damRegCoef_c; m_damRegCoef_c = NULL;
    }

    if (m_X1TimeRegCoef_c != NULL)
    {
      delete m_X1TimeRegCoef_c; m_X1TimeRegCoef_c = NULL;
    }

}


InverseVerticallyIntegratedVelocitySolver::InverseVerticallyIntegratedVelocitySolver()
  : m_config(), m_time(0.0), m_prev_time(-123.0)
{
  
}

InverseVerticallyIntegratedVelocitySolver::~InverseVerticallyIntegratedVelocitySolver()
{

  free(m_C);
  free(m_Cmasked);
  free(m_COrigin);
  free(m_bestC);
  free(m_muCoef);
  free(m_muCoefOrigin);
  free(m_damRegCoef);
  free(m_oneMinusMuSqBeta);
  free(m_X1TimeRegCoef);
  free(m_scaledX1Norm);
  free(m_lapC);
  free(m_gradCSq);
  free(m_lapMuCoef);
  free(m_gradMuCoefSq);
  free(m_lapX);
  free(m_gradXSq);
  free(m_velObs);
  free(m_velCoef);
  free(m_velb);
  free(m_velocityMisfit);
}

void
InverseVerticallyIntegratedVelocitySolver::Configuration::parse(const char* a_prefix)
{
  
  ParmParse pp(a_prefix);

  m_minLevelForOptimization = -1;
  pp.query("minLevelForOptimization", m_minLevelForOptimization);

  m_minTimeBetweenOptimizations = -1.0;
  pp.query("minTimeBetweenOptimizations", m_minTimeBetweenOptimizations  );

  m_dtTypical = 1.0/12.0;
  pp.query("dtTypical",m_dtTypical);
  
  m_CGmaxIter = 16;
  pp.query("CGmaxIter",m_CGmaxIter);
  
  m_CGtol = 1.0e-3;
  pp.query("CGtol",m_CGtol);

  m_CGhang = 0.999;
  pp.query("CGhang",m_CGhang);
  
  m_CGlineSearchMethod = CGOPT_LINE_SEARCH_SECANT;
  m_CGlineSearchInitialStep = 1.0e-7;
  m_CGlineSearchStepMaxGrow = 2.0;
  m_CGlineSearchMaxIter = 20;
  m_CGlineSearchTol = 1.0e-1;
  
  // backward compat.
  pp.query("CGsecantParameter",m_CGlineSearchInitialStep);
  pp.query("CGsecantStepMaxGrow",m_CGlineSearchStepMaxGrow);
  pp.query("CGsecantMaxIter",m_CGlineSearchMaxIter);
  pp.query("CGsecantTol",m_CGlineSearchTol);
  // new names
  pp.query("CGlineSearchMethod",m_CGlineSearchMethod);
  pp.query("CGlineSearchInitialStep",m_CGlineSearchInitialStep);
  pp.query("CGlineSearchStepMaxGrow",m_CGlineSearchStepMaxGrow);
  pp.query("CGlineSearchMaxIter",m_CGlineSearchMaxIter);
  pp.query("CGlineSearchTol",m_CGlineSearchTol);
  

  m_CGrestartInterval = 9999;
  pp.query("restartInterval",m_CGrestartInterval);

  m_velMisfitCoefficient = 1.0;
  pp.query("velMisfitCoefficient",m_velMisfitCoefficient);
 
  m_thicknessThreshold = 100.0;
  pp.query("thicknessThreshold",m_thicknessThreshold);
 
  {
    std::string s = "speed";
    pp.query("velMisfitType",s);
    if (s == "speed")
      {
	m_velMisfitType = speed;
      }
    else if (s == "velocity")
      {
	m_velMisfitType = velocity;
      }
    else if (s == "log_speed")
      {
	m_velMisfitType = log_speed;
      }
    else
      {
	MayDay::Error("unknown control.velMisfitType");
      }
  }
  
  // specify observed data
  m_velObs_c = SurfaceData::parse("control.velCoef");
  CH_assert(m_velObs_c != NULL);
  m_velObs_x = SurfaceData::parse("control.xVel");
  CH_assert(m_velObs_x != NULL);
  m_damRegCoef_c = SurfaceData::parse("control.damRegCoef");
  m_X1TimeRegCoef_c = SurfaceData::parse("control.X1TimeRegCoef");
#if CH_SPACEDIM > 1
  m_velObs_y = SurfaceData::parse("control.yVel");
  if ((m_velObs_y == NULL) && (m_velMisfitType == speed))
    {
      m_velObs_y = new ZeroData;
    }
#endif
  m_gradientFactor = SurfaceData::parse("control.gradientFactor");
  if (m_gradientFactor == NULL)
    {
      m_gradientFactor = new ConstantData(1.0);
    }

  // default: attempt to optimize w.r.t both X0 (C) and X1 (muCoef)
  m_optimizeX0 = true;
  pp.query("optimizeX0", m_optimizeX0 );
  m_optimizeX1 = true;
  pp.query("optimizeX1", m_optimizeX1 );
  CH_assert(m_optimizeX0 || m_optimizeX1); // at least one of these should be true

  m_fracInformedX1TimeReg = false;
  pp.query("fracInformedX1TimeReg", m_fracInformedX1TimeReg);

  // evaluate *unregluarized* gradient of J with respect to mucoef (or rather, X1) in the shelf only?
  m_gradMuCoefShelfOnly = false;
  pp.query("gradMuCoefShelfOnly", m_gradMuCoefShelfOnly); 

  m_gradCsqRegularization = 0.0;
  pp.query("gradCsqRegularization",m_gradCsqRegularization);

  m_gradMuCoefsqRegularization = m_gradCsqRegularization;
  pp.query("gradMuCoefsqRegularization",m_gradMuCoefsqRegularization);

  m_oneMinusMuSqRegularization = 0.0;
  pp.query("oneMinusMuSqRegularization",m_oneMinusMuSqRegularization);
  pout() << "oneMinusMuSqRegularization:" << m_oneMinusMuSqRegularization << std::endl;

  m_l_curve_param_filename = "l_curve_params.csv";
  pp.query("l_curve_param_filename",m_l_curve_param_filename);

  m_muRegularizationType = 1;
  pp.query("muRegularizationType", m_muRegularizationType );

  m_X0Regularization = 0.0;
  pp.query("X0Regularization",m_X0Regularization);

  m_X1Regularization = 0.0;
  pp.query("X1Regularization",m_X1Regularization);
  
  m_regularizeCGradient = 0;
  pp.query("regularizeCGradient", m_regularizeCGradient);

  m_X0TimeRegularization = 0.0;
  pp.query("X0TimeRegularization",m_X0TimeRegularization);

  m_X1TimeRegularization = 0.0;
  pp.query("X1TimeRegularization",m_X1TimeRegularization);

  m_gradX0sqRegularization = 0.0;
  pp.query("gradX0sqRegularization",  m_gradX0sqRegularization);

  m_gradX1sqRegularization = 0.0;
  pp.query("gradX1sqRegularization", m_gradX1sqRegularization);

  m_initialLowerC = 1.0e+0;
  pp.query("initialLowerC",m_initialLowerC);
  m_initialUpperC = 1.0e+6;
  pp.query("initialUpperC",m_initialUpperC);
  CH_assert(m_initialUpperC > m_initialLowerC);

  m_initialLowerMuCoef = 0.01;
  pp.query("initialLowerMuCoef",m_initialLowerMuCoef);
  m_initialUpperMuCoef = 100.0;
  pp.query("initialUpperMuCoef",m_initialUpperMuCoef);
  CH_assert(m_initialUpperMuCoef > m_initialLowerMuCoef);

  

  m_lowerX0 = - 3.0;
  m_upperX0 = + 3.0;
  pp.query("lowerX0",m_lowerX0);
  pp.query("upperX0",m_upperX0);
  CH_assert(m_upperX0 > m_lowerX0);


  m_lowerX1 = - 2.0;
  m_upperX1 = + 2.0;
  pp.query("lowerX1",m_lowerX1);
  pp.query("upperX1",m_upperX1);
  CH_assert(m_upperX1 > m_lowerX1);

  {
    std::string s = "none";
    pp.query("boundMethod",s);
    if ( (s == "none") || (s == "None"))
      {
	m_boundMethod = none;
      }
    else if ( (s == "projection") || (s == "Projection") )
      {
	m_boundMethod = projection;
      }
    else
      {
	MayDay::Error("unknown control.boundMethod");
      }
  }

  m_rankRedBlack = false;
  pp.query("rankRedBlack",m_rankRedBlack);
  

  m_writeInnerSteps = false;
  pp.query("writeInnerSteps",m_writeInnerSteps);

  m_innerStepFileNameBase = "ControlInner";
  pp.query("innerStepFileNameBase",m_innerStepFileNameBase);
  
  m_outerStepFileNameBase = "ControlOuter";
  pp.query("outerStepFileNameBase",m_outerStepFileNameBase);
  
  m_outerStepsNum = 0;
  pp.query("outerStepsNum",m_outerStepsNum);
  m_writeSelectOuterSteps = false;
  if (m_outerStepsNum > 0)
    {
      m_writeSelectOuterSteps = true;
      pp.getarr("outerSteps", m_outerSteps, 0, m_outerStepsNum);
    }
  
  // CG save and read options  
  m_CGsaveStateInterval = 987654321; // I guess this is more than anyone will do
  pp.query("CGsaveStateInterval", m_CGsaveStateInterval);
  m_CGsaveStateFileNameBase = "ControlCGState";
  pp.query("CGsaveStateFileNameBase",m_CGsaveStateFileNameBase);
  m_CGreadStateFile = "";
  pp.query("CGreadStateFile", m_CGreadStateFile);
  m_CGreadStateIter = -1;
  m_CGfreshStart = false;
  if (m_CGreadStateFile != "")
    {
      pp.get("CGreadStateIter",m_CGreadStateIter);
      pp.query("CGfreshStart",m_CGfreshStart);
    }
  
}


void 
InverseVerticallyIntegratedVelocitySolver::define
(const AmrIceBase& a_amrIce,
 const ProblemDomain& a_coarseDomain,
 ConstitutiveRelation* a_constitutiveRelation,
 BasalFrictionRelation* a_basalFrictionRelation,
 const Vector<DisjointBoxLayout>& a_grids,
 const Vector<int>& a_refRatio,
 const RealVect& a_dxCrse,
 IceThicknessIBC* a_thicknessIBC,
 int a_numLevels)
{
  m_optimization_done = false;
  m_basalFrictionRelation = a_basalFrictionRelation;
  m_constitutiveRelation = a_constitutiveRelation;
  m_thicknessIBC = a_thicknessIBC;
  m_amrIce = &a_amrIce;
  m_finest_level = a_numLevels - 1;
  m_grids = a_grids;
  m_vectOps.resize(a_numLevels);
  m_refRatio = a_refRatio,
    m_dx.resize(a_numLevels);
  m_dx[0] = a_dxCrse; 
  for (int lev = 1; lev < a_numLevels; lev++)
    {
      m_dx[lev] =  m_dx[lev-1] / Real(m_refRatio[lev-1]); 
    }
  m_domain.resize(a_numLevels);
  for (int lev = 0; lev < a_numLevels; lev++)
    {
      m_domain[lev] = m_grids[lev].physDomain();
    }

  m_velocityInitialised = false;
 
  // read configuration from ParmParse table
  m_config.parse();

  create(m_C,1,IntVect::Unit);
  create(m_Cmasked,1,IntVect::Unit);
  create(m_COrigin,1,IntVect::Unit);
  create(m_bestC,1,IntVect::Unit);

  create(m_muCoef,1,IntVect::Unit);
  create(m_muCoefOrigin,1,IntVect::Unit);
  create(m_bestMuCoef,1,IntVect::Unit);

  create(m_lapC,1,IntVect::Zero);
  create(m_gradCSq,1,IntVect::Zero);
 
  create(m_lapMuCoef,1,IntVect::Zero);
  create(m_gradMuCoefSq,1,IntVect::Zero);

  create(m_lapX,NXCOMP,IntVect::Zero);
  create(m_gradXSq,NXCOMP,IntVect::Zero);

  create(m_velObs,SpaceDim,IntVect::Unit);
  create(m_velCoef,1,IntVect::Unit);
  create(m_velocityMisfit,1,IntVect::Unit);

  create(m_damRegCoef,1,IntVect::Unit);
  create(m_oneMinusMuSqBeta,1,IntVect::Zero);

  create(m_X1TimeRegCoef,1,IntVect::Unit);
  create(m_scaledX1Norm,1,IntVect::Zero);

  create(m_velb,SpaceDim,IntVect::Unit);
  create(m_adjVel,SpaceDim,IntVect::Unit);
  create(m_adjRhs,SpaceDim,IntVect::Unit);

}

int InverseVerticallyIntegratedVelocitySolver::solve
(Vector<LevelData<FArrayBox>* >& a_horizontalVel,
 Vector<LevelData<FArrayBox>* >& a_calvedIce,
 Real& a_initialResidualNorm, Real& a_finalResidualNorm,
 const Real a_convergenceMetric,
 const Vector<LevelData<FArrayBox>* >& a_rhs,
 const Vector<LevelData<FArrayBox>* >& a_C,
 const Vector<LevelData<FArrayBox>* >& a_C0,
 const Vector<LevelData<FArrayBox>* >& a_A,
 const Vector<LevelData<FArrayBox>* >& a_muCoef,
 Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
 Real a_time,
 int a_lbase, int a_maxLevel)
{
  CH_TIME("InverseVerticallyIntegratedVelocitySolver::solve");
  pout() << "InverseVerticallyIntegratedVelocitySolver::solve "<< std::endl;

  m_bestMisfit = 1.23456789e+300;
  m_coordSys = a_coordSys;
  m_A = a_A;
  m_C0 = a_C0;
  m_rhs = a_rhs;
  m_calvedIce = a_calvedIce;
  m_time = a_time;

  bool skipOptimization = ((m_time - m_prev_time) < m_config.m_minTimeBetweenOptimizations);
  if (skipOptimization)
    {
      pout() << "InverseVerticallyIntegratedVelocitySolver::solve: skipping optimization" << std::endl;
    }

  
  //best fit velocity is to be output
  m_bestVel = a_horizontalVel;

  // assign input C to C_0 and limit if needed
  assign(m_COrigin, a_C);
  if (!skipOptimization)
    {
      for (int  lev = 0; lev < m_COrigin.size(); lev++)
	{
	  for (DataIterator dit(m_grids[lev]); dit.ok(); ++dit)
	    {
	      FArrayBox& c = (*m_COrigin[lev])[dit];
	      const Real& cmax = m_config.m_initialUpperC;
	      const Real& cmin = m_config.m_initialLowerC;
	      FORT_BOUNDCTRL(CHF_FRA1(c,0),
			     CHF_CONST_REAL(cmin),
			     CHF_CONST_REAL(cmax),
			     CHF_BOX(m_grids[lev][dit]));
	    } 
	}
    }
  
  if ( Abs(a_time - m_prev_time) < TINY_NORM) 
    {
      // if the time has not change since the last solve, 
      // this must be a refinement / regrid. Refinements
      // seem to work better if muCoef is not started
      // from the coarse mesh version (though maybe we
      // just need to interpolate differently?)
      pout() <<  " Optimization: muCoef <- 1 ";
      setToZero(m_muCoefOrigin);
      plus(m_muCoefOrigin,1.0);
    }
  else
    {
      // assign input muCoef to muCoef_0 and limit
      assign(m_muCoefOrigin, a_muCoef);
      if (!skipOptimization)
	{
	  for (int  lev = 0; lev < m_muCoefOrigin.size(); lev++)
	    {
	      for (DataIterator dit(m_grids[lev]); dit.ok(); ++dit)
		{
		  FArrayBox& mc = (*m_muCoefOrigin[lev])[dit];
		  const Real& mcmax = m_config.m_initialUpperMuCoef;
		  const Real& mcmin = m_config.m_initialLowerMuCoef;
		  FORT_BOUNDCTRL(CHF_FRA1(mc,0),
				 CHF_CONST_REAL(mcmin),
				 CHF_CONST_REAL(mcmax),
				 CHF_BOX(m_grids[lev][dit]));
		} 
	    }
	}
    }


  assign(m_velb, a_horizontalVel);
  
  Vector<LevelData<FArrayBox>* > X;
  create(X,2,IntVect::Unit);
  setToZero(X);

  //update observations
  for (int lev = 0; lev < m_velObs.size(); lev++)
    {
      LevelData<FArrayBox> u;
      // DFM (1/4/21) -- use D_TERM here to enable 1D (flowline) build
      D_TERM(
	     if (m_config.m_velObs_x)
	       {
		 aliasLevelData(u, m_velObs[lev], Interval(0,0) );
		 m_config.m_velObs_x->evaluate(u, *m_amrIce, lev, 0.0);
	       },
	     if (m_config.m_velObs_y)
	       {
		 aliasLevelData(u, m_velObs[lev], Interval(1,1) );
		 m_config.m_velObs_y->evaluate(u, *m_amrIce, lev, 0.0);
	       },
	     // D_TERM needs a 3D entry, which is nothing here.
             )

	if (m_config.m_velObs_c)
	  {
	    m_config.m_velObs_c->evaluate( *m_velCoef[lev], *m_amrIce, lev, 0.0);
	  }
      
      
      if (m_config.m_damRegCoef_c)
	{
	  m_config.m_damRegCoef_c->evaluate( *m_damRegCoef[lev], *m_amrIce, lev, 0.0);
	}
      
      if (m_config.m_X1TimeRegCoef_c)
	{
	  m_config.m_X1TimeRegCoef_c->evaluate( *m_X1TimeRegCoef[lev], *m_amrIce, lev, 0.0);
	}
      
    }
  
  if (X.size()  < m_config.m_minLevelForOptimization + 1)
    {
      // Avoid wasting time by not optimizing until the mesh has been
      // refined to some desired level. Instead, just fill C and phi
      // with the initial guess, and return the observed velocity
      // as the solution (which can then be used to refine as it was
      // in the older AMRIceControl class)
      mapX(X);
      assign(m_velb,m_velObs);
      assign(m_bestVel,m_velb);
      //      assign(m_bestC,m_Cmasked);
      assign(m_bestC,m_C);      
      assign(m_bestMuCoef, m_muCoef);
      m_optimization_done = false;
    }
  else
    {
      
      m_outerCounter = 0;
      m_innerCounter = 0;
      
      int CGmaxIter = m_config.m_CGmaxIter;
      //int CGminIter = m_config.m_CGminIter;
      if (skipOptimization)
	{
	  // just initialize the optimization, which means computing the first objective etc.
	  CGmaxIter = 0;
	  //CGminIter = -1;		
	}
  
      pout() << "Optimization: CGmaxIter = " << CGmaxIter << "  m_time = " << m_time << "  m_prev_time = " << m_prev_time  << std::endl;
      
      // attempt the optimization
      CGOptimize(*this,  X, CGmaxIter,
		 m_config.m_CGtol,
		 m_config.m_CGhang,
		 m_config.m_CGlineSearchMethod,
		 m_config.m_CGlineSearchInitialStep,
		 m_config.m_CGlineSearchStepMaxGrow, 
		 m_config.m_CGlineSearchMaxIter ,
		 m_config.m_CGlineSearchTol,
		 m_outerCounter);
      
      m_optimization_done = true;

      if (CGmaxIter > 0)
	m_prev_time = m_time;

    }

     
  free(X);

  return 0;
}

// duplicate storage of a_b in a_a
void InverseVerticallyIntegratedVelocitySolver::create
(Vector<LevelData<FArrayBox>* >& a_a, 
 const  Vector<LevelData<FArrayBox>* >& a_b)
{
  a_a.resize(a_b.size());
  for (int lev = 0; lev < a_b.size(); lev++)
    {
      if (a_a[lev] == NULL)
	a_a[lev] = new LevelData<FArrayBox>();
      m_vectOps[lev].create(*a_a[lev],*a_b[lev]);
    }
} 

//apply preconditioner s = M^{-1}r
void InverseVerticallyIntegratedVelocitySolver::preCond
(Vector<LevelData<FArrayBox>* >& a_s, 
 const Vector<LevelData<FArrayBox>* >& a_r)
{
  assign(a_s,a_r);
}



// set a_x = s * a_x
void InverseVerticallyIntegratedVelocitySolver::setToZero
(Vector<LevelData<FArrayBox>* >& a_x)
{

  for (int lev = 0; lev < a_x.size(); lev++)
    m_vectOps[lev].setToZero(*a_x[lev]);
}


// set a_x = a_x + a_r
void InverseVerticallyIntegratedVelocitySolver::plus
(Vector<LevelData<FArrayBox>* >& a_x, Real a_r)
{

  for (int lev = 0; lev < a_x.size(); lev++)
    m_vectOps[lev].plus(*a_x[lev], a_r);
}


// set a_x = s * a_x
void InverseVerticallyIntegratedVelocitySolver::scale
(Vector<LevelData<FArrayBox>* >& a_x, 
 const  Real a_s)
{

  for (int lev = 0; lev < a_x.size(); lev++)
    m_vectOps[lev].scale(*a_x[lev],a_s);
}

// set a_y = a_x
void InverseVerticallyIntegratedVelocitySolver::assign
(Vector<LevelData<FArrayBox>* >& a_y, 
 const Vector<LevelData<FArrayBox>* >& a_x)
{
  for (int lev = 0; lev < std::min(a_x.size(),a_y.size()); lev++)
    m_vectOps[lev].assign(*a_y[lev],*a_x[lev]);
}
  
// set a_y = a_y + a_s * a_x
void InverseVerticallyIntegratedVelocitySolver::incr
(Vector<LevelData<FArrayBox>* >& a_y, 
 const Vector<LevelData<FArrayBox>* >& a_x, Real a_s)
{
  for (int lev = 0; lev < std::min(a_x.size(),a_y.size()); lev++)
    m_vectOps[lev].incr(*a_y[lev],*a_x[lev],a_s);
}

// set a_y = a_y + a_s * a_x subject to bounds
void InverseVerticallyIntegratedVelocitySolver::incrState
(Vector<LevelData<FArrayBox>* >& a_y, 
 const Vector<LevelData<FArrayBox>* >& a_x, Real a_s)
{
  incr(a_y, a_x, a_s);
  for (int lev = 0; lev < std::min(a_x.size(),a_y.size()); lev++)
    {
      for (DataIterator dit(m_grids[lev]);dit.ok();++dit)
       {
  	 FArrayBox& y = (*a_y[lev])[dit];
	 
  	 FORT_BOUNDCTRL(CHF_FRA1(y,0),
  			CHF_CONST_REAL(m_config.m_lowerX0),
  			CHF_CONST_REAL(m_config.m_upperX0),
  			CHF_BOX(m_grids[lev][dit]));

  	 FORT_BOUNDCTRL(CHF_FRA1(y,1),
  			CHF_CONST_REAL(m_config.m_lowerX1),
  			CHF_CONST_REAL(m_config.m_upperX1),
  			CHF_BOX(m_grids[lev][dit]));
	 
       } 
    }
}
  


// return a_y.a_x
Real InverseVerticallyIntegratedVelocitySolver::dotProduct
(Vector<LevelData<FArrayBox>* >& a_1, 
 const Vector<LevelData<FArrayBox>* >& a_2)
{
  
  //AMR calculation copied and pasted from MultiLevelLinearOp

  CH_TIME("InverseVerticallyIntegratedVelocitySolver::dotProduct");

  // want to do this in an AMR way, so need to generate a temporary
  Vector<LevelData<FArrayBox>* > temp1, temp2;
  create(temp1, a_1);
  create (temp2, a_2);

  // first set to zero, then call assign, since assign only sets
  // valid regions (that way ghost cells are set to zero)

  setToZero(temp1);
  setToZero(temp2);

  assign(temp1, a_1);
  assign(temp2, a_2);

  // now set covered regions to zero
  for (int level =0 ; level<temp1.size()-1; level++)
    {
      LevelData<FArrayBox>& temp1Level = *temp1[level];
      LevelData<FArrayBox>& temp2Level = *temp2[level];

      CH_assert(temp1[level]->getBoxes() == temp2[level]->getBoxes());
      CH_assert(temp1[level+1] != NULL);

      int nRefFine = m_refRatio[level];
      const DisjointBoxLayout& finerGrids = temp1[level+1]->getBoxes();
      const DisjointBoxLayout& levelGrids = temp1[level]->getBoxes();
      DataIterator levelDit = levelGrids.dataIterator();
      LayoutIterator finerLit = finerGrids.layoutIterator();

      for (levelDit.begin(); levelDit.ok(); ++levelDit)
        {
          const Box& thisBox = levelGrids[levelDit];
          for (finerLit.begin(); finerLit.ok(); ++finerLit)
            {
              Box testBox = finerGrids[finerLit];
              testBox.coarsen(nRefFine);
              testBox &= thisBox;
              if (!testBox.isEmpty())
                {
                  temp1Level[levelDit].setVal(0.0, testBox,
                                              0, temp1Level.nComp());

                  temp2Level[levelDit].setVal(0.0, testBox,
                                              0, temp2Level.nComp());
                }
            } // end loop over finer boxes
        } // end loop over boxes on this level
    } // end loop over levels for setting covered regions to zero;

  // now loop over levels and call AMRLevelOp dotProduct
  Real prod = 0.0;
  for (int level = 0; level<temp1.size(); level++)
    {
      Real levelProd = m_vectOps[level].dotProduct(*temp1[level],
						   *temp2[level]);

      // incorporate scaling for AMR dot products
      RealVect dxLev = m_dx[level];
      Real scale = D_TERM(dxLev[0], *dxLev[1], *dxLev[2]);
      levelProd *= scale;
      prod += levelProd;
    }

  free(temp1);
  free(temp2);

  return prod;

}

// max no of CG iterations before restart.
int InverseVerticallyIntegratedVelocitySolver::nDoF(const Vector<LevelData<FArrayBox>* >& x)
{
  return m_config.m_CGrestartInterval;
}


/// convert a_x -> C, muCoef, and update dependent fields (e.g laplacians)
void 
InverseVerticallyIntegratedVelocitySolver::mapX(const Vector<LevelData<FArrayBox>* >& a_x)
{

  CH_TIME("InverseVerticallyIntegratedVelocitySolver::mapX");
  pout() << "InverseVerticallyIntegratedVelocitySolver::mapX" << std::endl;
  
  
  updateInvalid(a_x);

  if (m_config.m_rankRedBlack)
    {
      rankRedBlackSolution(a_x);
    }
  
  // convert a_x -> C, muCoef
  for (int lev=0; lev <= m_finest_level;lev++)
    {
      LevelData<FArrayBox>& levelX =  *a_x[lev];
      LevelData<FArrayBox>& levelC =  *m_C[lev];
      LevelData<FArrayBox>& levelCmasked =  *m_Cmasked[lev];
      LevelData<FArrayBox>& levelCOrigin =  *m_COrigin[lev];
      LevelData<FArrayBox>& levelMuCoef =  *m_muCoef[lev];
      LevelData<FArrayBox>& levelMuCoefOrigin =  *m_muCoefOrigin[lev];
      const DisjointBoxLayout levelGrids =  m_grids[lev];
      for (DataIterator dit(levelGrids);dit.ok();++dit)
	{
	  const Box& box = levelGrids[dit];
	  
	  levelC[dit].copy(levelCOrigin[dit]);
	  FORT_BOUNDEXPCTRL(CHF_FRA1(levelC[dit],0),
	  		    CHF_CONST_FRA1(levelX[dit],CCOMP),
	  		    CHF_CONST_REAL(m_config.m_lowerX0),
			    CHF_CONST_REAL(m_config.m_upperX0),
	  		    CHF_BOX(levelGrids[dit]));
	  
	  levelCmasked[dit].copy(levelC[dit]);

	  levelMuCoef[dit].copy(levelMuCoefOrigin[dit]);
	  FORT_BOUNDEXPCTRL(CHF_FRA1(levelMuCoef[dit],0),
	  		    CHF_CONST_FRA1(levelX[dit],MUCOMP),
	  		    CHF_CONST_REAL(m_config.m_lowerX1),
			    CHF_CONST_REAL(m_config.m_upperX1),
	  		    CHF_BOX(levelGrids[dit]));
	}
    }

  updateInvalid(m_C);
  updateInvalid(m_muCoef);

  for (int lev=0; lev <= m_finest_level;lev++)
    {

      LevelData<FArrayBox>& levelX =  *a_x[lev];
      LevelData<FArrayBox>& levelC =  *m_C[lev];
      LevelData<FArrayBox>& levelMuCoef =  *m_muCoef[lev];

      IceUtility::applyHelmOp(*m_lapC[lev], levelC, 0.0, 1.0,  m_grids[lev], m_dx[lev]);
      IceUtility::applyGradSq(*m_gradCSq[lev], levelC,m_grids[lev], m_dx[lev]);
      
      IceUtility::applyHelmOp(*m_lapMuCoef[lev], levelMuCoef, 0.0, 1.0, m_grids[lev], m_dx[lev]);
      IceUtility::applyGradSq(*m_gradMuCoefSq[lev], levelMuCoef,  m_grids[lev], m_dx[lev]);
      
      IceUtility::applyHelmOp(*m_lapX[lev], levelX, 0.0, 1.0,  m_grids[lev], m_dx[lev]);
      IceUtility::applyGradSq(*m_gradXSq[lev], levelX,m_grids[lev], m_dx[lev]);
    }
    

  {
    Real maxC = computeMax(m_C,m_refRatio);
    Real minC = computeMin(m_C,m_refRatio);
    Real maxMuCoef = computeMax(m_muCoef,m_refRatio);
    Real minMuCoef = computeMin(m_muCoef,m_refRatio);
    Real maxX0 = computeMax(a_x,m_refRatio,Interval(0,0));
    Real minX0 = computeMin(a_x,m_refRatio,Interval(0,0));
    Real maxX1 = computeMax(a_x,m_refRatio,Interval(1,1));
    Real minX1 = computeMin(a_x,m_refRatio,Interval(1,1));
    
    pout() 
      << "InverseVerticallyIntegratedVelocitySolver::mapX "
      << " max X[0] = " << maxX0
      << " min X[0] = " << minX0
      << " max X[1] = " << maxX1
      << " min X[1] = " << minX1 
      << " max C = " << maxC 
      << " min C = " << minC
      << " max muCoef = " << maxMuCoef
      << " min muCoef = " << minMuCoef
      << std::endl;
  }

  //Usual adjustments to C : set to 0 in shelves. add wall drag

  for (int lev=0; lev <= m_finest_level;lev++)
    {	
      IceUtility::setFloatingBasalFriction(*m_Cmasked[lev], *m_coordSys[lev], m_grids[lev]);
    }

  // add drag due to ice in contact with ice-free rocky walls
  ParmParse ppamr("amr");
  bool wallDrag = true; 
  ppamr.query("wallDrag", wallDrag);
  if (wallDrag)
    {
      Real wallDragExtra = 0.0;
      ppamr.query("wallDragExtra", wallDragExtra);
      
      for (int lev=0; lev <= m_finest_level;lev++)
	{
	  const LevelSigmaCS& levelCS = *m_coordSys[lev];
	  LevelData<FArrayBox>& levelCmasked = *m_Cmasked[lev];
	  const LevelData<FArrayBox>& levelC = *m_C[lev]; //not set to zero in shelves
	  const DisjointBoxLayout levelGrids = m_grids[lev];
	  for (DataIterator dit(levelGrids);dit.ok();++dit)
	    {
	      FArrayBox wallC(levelGrids[dit],1);
	      wallC.setVal(0.0);
	      IceUtility::addWallDrag(wallC, 
				      levelCS.getFloatingMask()[dit], levelCS.getSurfaceHeight()[dit],
				      levelCS.getH()[dit], levelCS.getTopography()[dit], 
				      levelC[dit], wallDragExtra,m_dx[lev],levelGrids[dit]);
	      
	      levelCmasked[dit] += wallC;
	    }
	}
    }
}

void 
InverseVerticallyIntegratedVelocitySolver::solveStressEqn
(Vector<LevelData<FArrayBox>* >& a_u,
 const bool a_adjoint,
 const Vector<LevelData<FArrayBox>* >& a_rhs,
 const Vector<LevelData<FArrayBox>* >& a_C,
 const Vector<LevelData<FArrayBox>* >& a_C0,
 const Vector<LevelData<FArrayBox>* >& a_A,
 const Vector<LevelData<FArrayBox>* >& a_muCoef)
{

  CH_TIME("InverseVerticallyIntegratedVelocitySolver::solveStressEqn");
  pout() << "InverseVerticallyIntegratedVelocitySolver::solveStressEqn" << std::endl;
  
  
  JFNKSolver jfnkSolver;
  jfnkSolver.define(m_domain[0], m_constitutiveRelation , m_basalFrictionRelation,
		    m_grids, m_refRatio, m_dx[0], m_thicknessIBC, m_finest_level+1);
    
  Real initialNorm = 1.0; Real finalNorm = 1.0; Real convergenceMetric=-1.0;

  if (a_adjoint)
    {
      JFNKSolver::Configuration cfg 
	= jfnkSolver.get_evil_configuration();

      cfg.m_maxRelaxIter = cfg.m_maxRelaxIter * 2;
      cfg.m_RelaxHang = 0.25;
 
    }
  
  bool linear = a_adjoint;

  
  jfnkSolver.solve(a_u, m_calvedIce, 
		   initialNorm, finalNorm, convergenceMetric, linear, 
		   a_rhs, a_C, a_C0, a_A, a_muCoef, m_coordSys, 0.0 , 0, m_finest_level);


  if (a_adjoint)
    {
      pout() << " ... adjoint equation final residual = " << finalNorm  << "/" << initialNorm << " = " << finalNorm/initialNorm << std::endl;
      //CH_assert(finalNorm <= initialNorm);
    }
}

void 
InverseVerticallyIntegratedVelocitySolver::computeObjectiveAndGradient
(Real& a_fm, Real& a_fp, Vector<LevelData        <FArrayBox>* >& a_g, 
 const Vector<LevelData<FArrayBox>* >& a_x, int a_iter, bool a_inner)
{

  CH_TIME("InverseVerticallyIntegratedVelocitySolver::computeObjectiveAndGradient");
  pout() << "InverseVerticallyIntegratedVelocitySolver::computeObjectiveAndGradient" << std::endl;
  
  a_fm = 0.0;
  a_fp = 0.0;
  setToZero(a_g);

  //convert a_x -> C, muCoef
  mapX(a_x);

  //solve the forward problem (to update the effective viscosity & drag)
  if (a_inner)
    {
      // do not depend on the previous inner iteration
      assign(m_velb, m_bestVel);
    }
  solveStressEqn(m_velb,false,m_rhs,m_Cmasked,m_C0,m_A,m_muCoef);

  //todo : compute the L1L2 surface velocity ? was not too helpful before
  m_vels = m_velb;

  //construct rhs for the adjoint problem (and misfits)
  computeAdjointRhs();

  //compute objective function 
  Real vobj = computeSum(m_velocityMisfit, m_refRatio, m_dx[0][0]);
  Real sumGradCSq = computeSum(m_gradCSq,m_refRatio, m_dx[0][0]);
  Real sumGradMuSq = computeSum(m_gradMuCoefSq,m_refRatio, m_dx[0][0]);
  Real sumGradX0Sq =  computeSum(m_gradXSq,m_refRatio, m_dx[0][0], Interval(0,0));
  Real sumGradX1Sq =computeSum(m_gradXSq,m_refRatio, m_dx[0][0], Interval(1,1));
  Real normX0 = computeNorm(a_x,m_refRatio, m_dx[0][0], Interval(0,0));
  Real normX1 = computeNorm(a_x,m_refRatio, m_dx[0][0], Interval(1,1));
  
  Real sumOneMinusMuSqBeta = 0.0;
  if ( m_config.m_muRegularizationType == 2 )
    {
      pout() << "COMPUTING MU ONE MINUS MU SQUARED BETA" << std::endl;
      
      for (int lev = 0; lev <= m_finest_level; lev++)
	{
	  LevelData<FArrayBox>& oneMinusMuSqBeta = (*m_oneMinusMuSqBeta[lev]);
	  for (DataIterator dit(m_grids[lev]); dit.ok(); ++dit)
	    {
	      
	      oneMinusMuSqBeta[dit].setVal(1.0);
	      
	      const FArrayBox& muCoef = (*m_muCoef[lev])[dit];
	      const FArrayBox& drc = (*m_damRegCoef[lev])[dit];
	      
	      for ( BoxIterator bit(m_grids[lev][dit]); bit.ok(); ++bit ){
		const IntVect& iv = bit();
		oneMinusMuSqBeta[dit](iv)-= muCoef(iv);
		oneMinusMuSqBeta[dit](iv)*= oneMinusMuSqBeta[dit](iv);
		oneMinusMuSqBeta[dit](iv)*= drc(iv);
	      }
	    }
	}
      sumOneMinusMuSqBeta = computeSum(m_oneMinusMuSqBeta, m_refRatio, m_dx[0][0]);
    }
  

  Real sumScaledX1Norm = normX1*normX1;
  if ( m_config.m_fracInformedX1TimeReg )
    {
      sumScaledX1Norm = 0.0;
      pout() << "ADDING SCALED X1 NORM TO PENALTY TERM" << std::endl;
      for (int lev = 0; lev <= m_finest_level; lev++)
	{
	  LevelData<FArrayBox>& scaledX1Norm = (*m_scaledX1Norm[lev]);
	  for (DataIterator dit(m_grids[lev]);dit.ok();++dit)
	    {
	      scaledX1Norm[dit].setVal(0.0);
	      const FArrayBox& X = (*a_x[lev])[dit];
	      const FArrayBox& X1trc = (*m_X1TimeRegCoef[lev])[dit];
	      
	      for ( BoxIterator bit(m_grids[lev][dit]); bit.ok(); ++bit ){
		const IntVect& iv = bit();
		scaledX1Norm[dit](iv)+= X1trc(iv);
		scaledX1Norm[dit](iv)*= X(iv,MUCOMP);
		scaledX1Norm[dit](iv)*= X(iv,MUCOMP);
	      }
	    }
	}
      sumScaledX1Norm += computeSum(m_scaledX1Norm, m_refRatio, m_dx[0][0]);
    }

  a_fm = vobj ;
  a_fp =  m_config.m_gradCsqRegularization * sumGradCSq
    + m_config.m_gradMuCoefsqRegularization * sumGradMuSq
    + m_config.m_gradX0sqRegularization * sumGradX0Sq
    + m_config.m_gradX1sqRegularization * sumGradX1Sq
    + X0Regularization() * normX0*normX0
    + X1Regularization() * sumScaledX1Norm;
    + m_config.m_oneMinusMuSqRegularization * sumOneMinusMuSqBeta;
 
  pout() << ((a_inner)?"inner: ":"outer: ")
	 << " ||velocity misfit||^2 = " << vobj 
	 << " || grad C ||^2 = " << sumGradCSq
	 << " || grad muCoef ||^2 = " << sumGradMuSq
	 << " || X0 ||^2 = " << normX0*normX0
	 << " || X1 ||^2 = " << sumScaledX1Norm
	 << " || grad X0 ||^2 = " << sumGradX0Sq
	 << " || grad X1 ||^2 = " << sumGradX1Sq ;
	if ( m_config.m_muRegularizationType == 2 ){
		pout() << " || 1-Phi ||^2 * beta(x,y) = " << sumOneMinusMuSqBeta;
	}
  pout() << std::endl;
  

  if ( a_fm < m_bestMisfit && (!a_inner) )
    {
      //save the velocity, muCoef, and C;
      assign(m_bestVel,m_velb);
      assign(m_bestC,m_C);      
      assign(m_bestMuCoef, m_muCoef);
      m_bestMisfit = a_fm;
    }

  //solve the adjoint problem
  pout() << " solving adjoint equations... " << std::endl;
  //adjoint stress balance is linear, but we need to start from m_velb
  //to get the correct effective viscosity and basal drag
  assign(m_adjVel,m_velb);
  //attempt to avoid occasional divergence in shelf
  Real adjReg = 1.0;
  plus(m_Cmasked,adjReg);
  // Scaling C by m_basalFrictionRelation->power() results in a better approximation to the adjoint stress balance
  // there are more efficient wasy to achive this, but this is the simplest
  scale(m_Cmasked, m_basalFrictionRelation->power() );
  solveStressEqn(m_adjVel,true,m_adjRhs,m_Cmasked,m_C0,m_A,m_muCoef);
  // Put C back the way it was
  scale(m_Cmasked, 1.0/m_basalFrictionRelation->power());
  plus(m_Cmasked,-adjReg);

  
  //compute gradient 
  setToZero(a_g);
  
  //compute directional derivatives of unregularized problrm
  computeGradient(a_g, a_x);
  
  //add Tikhonov regularization
  regularizeGradient(a_g, a_x);

 
  //red-black scheme modifies gradient
  if (m_config.m_rankRedBlack)
    {
      rankRedBlackGradient(a_g);
    }
  
  if (m_config.m_boundMethod == Configuration::projection)
    applyProjection(a_g, a_x);

  // apply non-uniform factor 
  for (int lev = 0; lev <= m_finest_level; lev++)
    {
      LevelData<FArrayBox> f (m_grids[lev], 1, IntVect::Zero);
      m_config.m_gradientFactor->evaluate(f, *m_amrIce, lev, 0.0);
      for (DataIterator dit(m_grids[lev]); dit.ok(); ++dit)
	{
	  FArrayBox& g = (*a_g[lev])[dit];
	  for (int n = 0; n < g.nComp(); n++)
	    {
	      g.mult(f[dit], 0, n, 1);
	    }  
	}
    }

  //dump data
  if (m_config.m_writeInnerSteps)
    {
      writeState(innerStateFile(), m_innerCounter, a_x, a_g);
      m_innerCounter++;
    }
  if (!a_inner)
    {
      m_outerCounter = a_iter;  
      if (!m_config.m_writeSelectOuterSteps || (std::find(m_config.m_outerSteps.begin(), m_config.m_outerSteps.end(), m_outerCounter) != m_config.m_outerSteps.end()))
	{
	  writeState(outerStateFile(), m_innerCounter, a_x, a_g);
	}
      //m_outerCounter++;
    }
  
}

/// compute adjoint eqn rhs, plus related quantities (e.g misfits)
void 
InverseVerticallyIntegratedVelocitySolver::computeAdjointRhs()
{
  
 CH_TIME("InverseVerticallyIntegratedVelocitySolver::computeAdjointRhs");
 pout() << " InverseVerticallyIntegratedVelocitySolver::computeAdjointRhs" << std::endl;
 
 
  // rhs contribution due to velocity mismatch
  for (int lev=0; lev <= m_finest_level ;lev++)
    {
      for (DataIterator dit(m_grids[lev]);dit.ok();++dit)
	{
	  FArrayBox& adjRhs = (*m_adjRhs[lev])[dit];
	  FArrayBox& misfit = (*m_velocityMisfit[lev])[dit];
	  const FArrayBox& um = (*m_vels[lev])[dit];
	  const FArrayBox& uo = (*m_velObs[lev])[dit];
	  const FArrayBox& h = m_coordSys[lev]->getH()[dit];
	  FArrayBox& uc = (*m_velCoef[lev])[dit];
	  const Box& box = m_grids[lev][dit];
	  
	  if (m_config.m_velMisfitType == Configuration::speed)
	    {
	      FORT_ADJRHSSPEEDCTRL(CHF_FRA1(adjRhs,0), CHF_FRA1(adjRhs,1),
				   CHF_FRA1(misfit,0),
				   CHF_CONST_FRA1(um,0), CHF_CONST_FRA1(um,1),
				   CHF_CONST_FRA1(uo,0), CHF_CONST_FRA1(uo,1),
				   CHF_BOX(box));
	    }
	  else if (m_config.m_velMisfitType == Configuration::velocity)
	    {
	      FORT_ADJRHSVELCTRL(CHF_FRA1(adjRhs,0), CHF_FRA1(adjRhs,1),
				 CHF_CONST_FRA1(misfit,0),
				 CHF_CONST_FRA1(um,0), CHF_CONST_FRA1(um,1),
				 CHF_CONST_FRA1(uo,0), CHF_CONST_FRA1(uo,1),
				 CHF_BOX(box));
	    }
	  if (m_config.m_velMisfitType == Configuration::log_speed)
	    {
	      FORT_ADJRHSLOGSPDCTRL(CHF_FRA1(adjRhs,0), CHF_FRA1(adjRhs,1),
				    CHF_FRA1(misfit,0),
				    CHF_CONST_FRA1(um,0), CHF_CONST_FRA1(um,1),
				    CHF_CONST_FRA1(uo,0), CHF_CONST_FRA1(uo,1),
				    CHF_BOX(box));
	    }
	  else
	    {
	      CH_assert(m_config.m_velMisfitType < Configuration::MAX_VELOCITY_MISFIT_TYPE);
	    }

	  for (BoxIterator bit(box);bit.ok();++bit)
	    {
	      const IntVect& iv = bit();
	      if (uc(iv) < 0.975) uc(iv) = 0.0; 
	      if (h(iv) < m_config.m_thicknessThreshold) uc(iv) = 0.0;
	    }

	  for (int dir = 0; dir < SpaceDim; dir++)
	    {
	      adjRhs.mult(uc,0,dir);
	    }
	  misfit *= uc;
	   
	  adjRhs *= m_config.m_velMisfitCoefficient;
	  misfit *= m_config.m_velMisfitCoefficient;

	}
    }
  // rhs contribution due to velocity mismatch
}
 

/// unregularized gradient
void InverseVerticallyIntegratedVelocitySolver::computeGradient
(Vector<LevelData<FArrayBox>* >& a_g, 
 const  Vector<LevelData<FArrayBox>* >& a_x)
{
  
  CH_TIME("InverseVerticallyIntegratedVelocitySolver::computeGradient");
  pout() << " InverseVerticallyIntegratedVelocitySolver::computeGradient" << std::endl;
  
  if (m_config.m_optimizeX0)
      {
	// grad w.r.t x_0 (basal friction)  = - adjVel * vel * C 
	for (int lev = 0; lev <= m_finest_level; lev++)
	  {
	    for (DataIterator dit(m_grids[lev]);dit.ok();++dit)
	      {
		
		FArrayBox& G = (*a_g[lev])[dit];
		const FArrayBox& C = (*m_Cmasked[lev])[dit];
		const FArrayBox& u = (*m_velb[lev])[dit];
		const FArrayBox& lambda = (*m_adjVel[lev])[dit];
		FArrayBox alpha(u.box(), 1);
		
		Real scale = 1.0 ; // check!
		m_basalFrictionRelation->computeAlpha(alpha, u, C, scale,
						      *m_coordSys[lev],
						      dit, lev, alpha.box());
		
		
		FArrayBox t(G.box(),1);
		FORT_CADOTBCTRL(CHF_FRA1(t,0),
				CHF_CONST_FRA1(alpha,0),
				CHF_CONST_FRA(u),
				CHF_CONST_FRA(lambda),
				CHF_BOX(t.box()));
		t *= -1.0;
		G.plus(t,0,CCOMP);
	      }
	    
	  }
      } // end if m_optimizeX1
  
  if (m_config.m_optimizeX1)
    {
      // grad w.r.t x_1 (mu coef)  
      Vector<LevelData<FluxBox>*> faceVT;
      create(faceVT, SpaceDim, IntVect::Unit);
      Vector<LevelData<FluxBox>*> faceA;
      create(faceA, m_A[0]->nComp(), IntVect::Unit);
      for (int lev = 0; lev <= m_finest_level; lev++)
	{
	  CellToEdge(*m_A[lev], *faceA[lev]);
	}
      
      IceNonlinearViscousTensor state(m_grids, m_refRatio, m_domain, m_dx,
				      m_coordSys, m_velb, m_Cmasked, m_C0, m_finest_level ,
				      *m_constitutiveRelation, *m_basalFrictionRelation,
				      *m_thicknessIBC, m_A, faceA, 0.0, 0.0, 0, 0.0);
      state.setState(m_velb);
      state.computeViscousTensorFace(faceVT);
      
      
      for (int lev = 0; lev <= m_finest_level; lev++)
	{
	  for (DataIterator dit(m_grids[lev]);dit.ok();++dit)
	    { 
	      FArrayBox& G = (*a_g[lev])[dit];
	      const Box& box = m_grids[lev][dit]; 
	      const FArrayBox& lambda = (*m_adjVel[lev])[dit];
	      const FArrayBox& muCoef = (*m_muCoef[lev])[dit];
	      const FArrayBox& h = m_coordSys[lev]->getH()[dit];
	      const FluxBox& vt = (*faceVT[lev])[dit];
	      //need to avoid the calving front, so construct an interior mask 
	      BaseFab<bool> interior(box,1);
	      const Real tol = 1.0;
	      for (BoxIterator bit(box);bit.ok();++bit)
		{
		  const IntVect& iv = bit();
		  interior(iv) = (h(iv) > tol);
		  for (int dir = 0; dir < SpaceDim; dir++)
		    {
		      interior(iv) &= ( h(iv+BASISV(dir)) > tol);
		      interior(iv) &= ( h(iv-BASISV(dir)) > tol);
		    }
		}
	      
	      FArrayBox t(box,1);
	      t.setVal(0.0);
	      for (int facedir = 0; facedir < SpaceDim; facedir++)
		{
		  const IntVect e = BASISV(facedir);
		  for (BoxIterator bit(box);bit.ok();++bit)
		    {
		      for (int dir = 0; dir < SpaceDim; dir++)
			{
			  const IntVect& iv = bit();
			  if ( interior(iv) )
			    {
			      t(iv) += (lambda(iv+e,dir)-lambda(iv,dir))* vt[facedir](iv+e,dir)
				+ (lambda(iv,dir)-lambda(iv-e,dir))* vt[facedir](iv,dir);
			    }
			}
		    }
		}
	      
	      RealVect oneOnTwoDx = 1.0/ (2.0 * m_dx[lev]);
	      t.mult(-oneOnTwoDx[0]);
	      t.mult(muCoef);

	      // restrict X1 gradient to shelf only?
	      if (m_config.m_gradMuCoefShelfOnly)
		{
		  const BaseFab<int>& mask = m_coordSys[lev]->getFloatingMask()[dit];
		  for (BoxIterator bit(box);bit.ok();++bit)
		    {
		      const IntVect& iv = bit();
		      if (mask(iv) != FLOATINGMASKVAL)
			{
			  t(iv) = 0.0;
			}
		    }
		}
	      
	      G.plus(t,0,MUCOMP);
	    }
	}
      free(faceVT);
      free(faceA);
    } // end if m_config.m_optimizeX1
}

/// update invalid cells (ghost + covered coarse cells)
void InverseVerticallyIntegratedVelocitySolver::updateInvalid(const Vector<LevelData<FArrayBox>* >& a_x)
{
  CH_TIME("InverseVerticallyIntegratedVelocitySolver::updateInvalid");
  pout() << " InverseVerticallyIntegratedVelocitySolver::updateInvalid" << std::endl;

  
  for (int lev=m_finest_level; lev > 0 ;lev--)
    {
      CoarseAverage avg(m_grids[lev],a_x[lev]->nComp(),m_refRatio[lev-1]);
      avg.averageToCoarse(*a_x[lev-1],*a_x[lev]);
    }
  
  for (int lev=1; lev <= m_finest_level;lev++)
    {
      // coarse-finery
      Real time_interp_coeff = 0.0;
      int nGhost = 1;
      int ncomp = (*a_x[lev]).nComp();
      PiecewiseLinearFillPatch lin(m_grids[lev],  m_grids[lev-1],ncomp,          
				   m_domain[lev-1],m_refRatio[lev-1], nGhost);
      lin.fillInterp(*a_x[lev], *a_x[lev-1],*a_x[lev-1],time_interp_coeff, 0, 0, ncomp);   
    }

  for (int lev=0; lev <= m_finest_level;lev++)
    {
      a_x[lev]->exchange();
      // reflection at the boundaries
      for (int dir = 0; dir < SpaceDim; dir++)
	{
	  if (! m_domain[lev].isPeriodic(dir))
	    {
	      ReflectGhostCells(*a_x[lev], m_domain[lev], dir, Side::Lo);
	      ReflectGhostCells(*a_x[lev], m_domain[lev], dir, Side::Hi);
	    }
	}
    }
}


/** modify the solution vector x e red-black scheme 
    x(i,j) <- 1/4 ( x(i+1,j) + x(i-1,j) + ... )/4 and g(i,j) = 0
**/
void InverseVerticallyIntegratedVelocitySolver::rankRedBlackSolution
(const Vector<LevelData<FArrayBox>* >& a_x)
{

  CH_TIME("InverseVerticallyIntegratedVelocitySolver::rankRedBlackSolution");
  pout() << " InverseVerticallyIntegratedVelocitySolver::rankRedBlackSolution" << std::endl;
  
  // rank reduction scheme: interpolate x0 from 'red' to 'black' and x1 from 'black' to 'red'
  for (int lev=0; lev <= m_finest_level;lev++)
    {
      const DisjointBoxLayout levelGrids =  m_grids[lev];
      for (DataIterator dit(levelGrids);dit.ok();++dit)
	{
	  const Box& box = levelGrids[dit];
	  FArrayBox& x = (*a_x[lev])[dit];
	  int comp;
	  // \todo support 1D, make efficient.
	  for (BoxIterator bit(box); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      IntVect e = iv + BASISV(0);
	      IntVect w = iv - BASISV(0);
	      IntVect s = iv - BASISV(1);
	      IntVect n = iv + BASISV(1);
	      int comp = (iv[0]%2 == iv[1]%2)?CCOMP:MUCOMP; // \todo fix 2D only
	      x(iv,comp) = 0.25 * (x(e,comp) + x(w,comp) + x(n,comp) + x(s,comp));
	    }
	}
    }
}

/** modify the gradient vector according to the red-black scheme 
    
    g(i,j) <- g(i,j) + ( g(i+1,j) + g(i-1,j) + ... )/4

    \todo support 1D, make efficient.
*/
void InverseVerticallyIntegratedVelocitySolver::rankRedBlackGradient
(const Vector<LevelData<FArrayBox>* >& a_g)
{

  CH_TIME("InverseVerticallyIntegratedVelocitySolver::rankRedBlackGradient");
  pout() << " InverseVerticallyIntegratedVelocitySolver::rankRedBlackGradient" << std::endl;
  
  for (int lev = 0; lev <= m_finest_level; lev++)
    {

      if (lev > 0)
	{
	  //coarse fine ghosts
	  Real time_interp_coeff = 0.0;
	  int nGhost = 1;
	  int nComp = a_g[lev]->nComp();
	  PiecewiseLinearFillPatch lin(m_grids[lev],  m_grids[lev-1],nComp,          
				       m_domain[lev-1],m_refRatio[lev-1], nGhost);
	  lin.fillInterp(*a_g[lev], *a_g[lev-1],*a_g[lev-1],time_interp_coeff, 0, 0, nComp);   
	}
      //fine-fine ghosts
      a_g[lev]->exchange();
      //boundary ghosts
      for (int dir = 0; dir < SpaceDim; dir++)
	{
	  if (! m_domain[lev].isPeriodic(dir))
	    {
	      ReflectGhostCells(*a_g[lev], m_domain[lev], dir, Side::Lo);
	      ReflectGhostCells(*a_g[lev], m_domain[lev], dir, Side::Hi);
	    }
	}
      
      const DisjointBoxLayout& grids = a_g[lev]->disjointBoxLayout(); 
      for (DataIterator dit(grids);dit.ok();++dit)
	{
	  FArrayBox& g = (*a_g[lev])[dit];
	  for (BoxIterator bit(grids[dit]); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      IntVect e = iv + BASISV(0);
	      IntVect w = iv - BASISV(0);
	      IntVect s = iv - BASISV(1);
	      IntVect n = iv + BASISV(1);
	      int comp = (iv[0]%2 == iv[1]%2)?MUCOMP:CCOMP;
	      g(iv,comp) += 0.25 * (g(e,comp) + g(w,comp) + g(n,comp) + g(s,comp));
	    }

	  for (BoxIterator bit(grids[dit]); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      IntVect e = iv + BASISV(0);
	      IntVect w = iv - BASISV(0);
	      IntVect s = iv - BASISV(1);
	      IntVect n = iv + BASISV(1);
	      int comp = (iv[0]%2 == iv[1]%2)?CCOMP:MUCOMP;
	      g(iv,comp) = 0.0; //needed? 
	    }
	  
	}

      int dbg=0;dbg++;
    }
}

/// add regularization terms R = (- a C lap(C)) etc to the gradient
void InverseVerticallyIntegratedVelocitySolver::regularizeGradient
(Vector<LevelData<FArrayBox>* >& a_g,  const  Vector<LevelData<FArrayBox>* >& a_x)
{

  CH_TIME("InverseVerticallyIntegratedVelocitySolver::regularizeGradient");
  pout() << " InverseVerticallyIntegratedVelocitySolver::regularizeGradient" << std::endl;

  
  for (int lev = 0; lev <= m_finest_level; lev++)
    {
      for (DataIterator dit(m_grids[lev]);dit.ok();++dit)
	{
	  FArrayBox& G = (*a_g[lev])[dit];
	  const FArrayBox& X = (*a_x[lev])[dit];
	  const FArrayBox& lapC = (*m_lapC[lev])[dit];
	  const FArrayBox& C = (*m_C[lev])[dit];
	  const FArrayBox& lapMuCoef = (*m_lapMuCoef[lev])[dit];
	  const FArrayBox& muCoef = (*m_muCoef[lev])[dit];
      const FArrayBox& drc = (*m_damRegCoef[lev])[dit];
      const FArrayBox& X1trc = (*m_X1TimeRegCoef[lev])[dit];
	  
	  FArrayBox t(m_grids[lev][dit],1);
	  // terms arising from (grad C)^2 penalty
	  if ( (m_config.m_gradCsqRegularization > 0.0) && m_config.m_optimizeX0)
	    {
	      t.copy(lapC);t*= C; t*= -m_config.m_gradCsqRegularization;
	      G.plus(t,0,CCOMP);
	    }

	  // terms arising from (grad muCoef)^2 penalty
	  if ( (m_config.m_gradMuCoefsqRegularization > 0.0) && m_config.m_optimizeX1)
	    {
	      t.copy(lapMuCoef);t*= muCoef; t*= -m_config.m_gradMuCoefsqRegularization;
	      G.plus(t,0,MUCOMP);
	    }
      // terms arising from damRegCoef*(1-\phi)^2 penalty
      if ( (m_config.m_muRegularizationType == 2) && (m_config.m_oneMinusMuSqRegularization > 0.0) && m_config.m_optimizeX1 )
        {
          FArrayBox minusOneMinusMuCoef(m_grids[lev][dit],1);
          minusOneMinusMuCoef.setVal(-1.0);
          minusOneMinusMuCoef+= muCoef;

          t.copy(muCoef); t*= minusOneMinusMuCoef; t*= drc; t*= m_config.m_oneMinusMuSqRegularization;
          G.plus(t,0,MUCOMP);

          if ( m_config.m_gradMuCoefsqRegularization > 0.0 )
            {
              // and adding the Tik term back in where we have fractures:
              FArrayBox invDRC(m_grids[lev][dit],1);
              invDRC.setVal(1.0);
              invDRC-=drc;

              if ( m_config.m_gradMuCoefsqRegularization > 0.0 ){
                  t.copy(lapMuCoef); t*= muCoef; t*= -m_config.m_gradMuCoefsqRegularization; t*=invDRC;
                  G.plus(t,0,MUCOMP);
              }
            }
        }

	  // terms arising from (X0)^2 penalty
	  if ( (X0Regularization() > 0.0) && m_config.m_optimizeX0)
	    {
	      t.copy(X,CCOMP,0);
	      t *= X0Regularization();
	      G.plus(t,0,CCOMP);
	    }

	  // terms arising from (X1)^2 penalty
      if ( (X1Regularization() > 0.0) && m_config.m_optimizeX1)
        {
          t.copy(X,MUCOMP,0);
          t *= X1Regularization();
          pout() << X1Regularization() << std::endl;
          if ( m_config.m_fracInformedX1TimeReg )
            {
              t *= X1trc;
            }
          G.plus(t,0,MUCOMP);
        }

	  // terms arising from (grad X0)^2 penalty
	  if ((m_config.m_gradX0sqRegularization > 0.0) && m_config.m_optimizeX0)
	    {
	      t.copy((*m_lapX[lev])[dit],CCOMP,0); t*= -m_config.m_gradX0sqRegularization;
	      CH_assert(t.norm() < 1.2345678e+300);
	      G.plus(t,0,CCOMP);
	    }

	  // terms arising from (grad X1)^2 penalty
	  if ((m_config.m_gradX1sqRegularization > 0.0) && m_config.m_optimizeX1)
	    {
	      t.copy((*m_lapX[lev])[dit],MUCOMP,0);
	      CH_assert(t.norm() < 1.2345678e+300);
	      t*= -m_config.m_gradX1sqRegularization;
	      G.plus(t,0,MUCOMP); 
	    }

	}
    }
}

/// set gradient a_g so that if a_x is at/outside bounds, descent direction (-a_g) is zero or inward
void InverseVerticallyIntegratedVelocitySolver::applyProjection
(Vector<LevelData<FArrayBox>* >& a_g, const  Vector<LevelData<FArrayBox>* >& a_x)
{

 CH_TIME("InverseVerticallyIntegratedVelocitySolver::applyProjection");
 pout() << "InverseVerticallyIntegratedVelocitySolver::applyProjection" << std::endl;
  
 for (int lev = 0; lev <= m_finest_level; lev++)
   {
     for (DataIterator dit(m_grids[lev]);dit.ok();++dit)
       {
	 FArrayBox& G = (*a_g[lev])[dit];
	 const FArrayBox& X = (*a_x[lev])[dit];
	 
	 FORT_HARDPOINTINCTRL(CHF_FRA1(G,CCOMP),
			      CHF_CONST_FRA1(X,CCOMP),
			      CHF_CONST_REAL(m_config.m_lowerX0),
			      CHF_CONST_REAL(m_config.m_upperX0),
			      CHF_BOX(G.box()));
	 
	 FORT_HARDPOINTINCTRL(CHF_FRA1(G,MUCOMP),
			      CHF_CONST_FRA1(X,MUCOMP),
			      CHF_CONST_REAL(m_config.m_lowerX1),
			      CHF_CONST_REAL(m_config.m_upperX1),
			      CHF_BOX(G.box()));
	 
	 
       }
   }
 
}


void InverseVerticallyIntegratedVelocitySolver::saveCGState
( int a_iter,
    const Vector<LevelData<FArrayBox>* >& a_x,
    const Vector<LevelData<FArrayBox>* >& a_r,
    const Vector<LevelData<FArrayBox>* >& a_s,
    const Vector<LevelData<FArrayBox>* >& a_d)
{

  if (a_iter % m_config.m_CGsaveStateInterval == 0)
    {
      std::string fileName(CGstateFile(a_iter));
	
      pout() << " InverseVerticallyIntegratedVelocitySolver::saveCGState: saving CG state to "
	     << fileName << std::endl;

      Vector<std::string> names(8+SpaceDim);
      int i(0);
      names[i++] = "x0"; names[i++] = "x1";
      names[i++] = "r0"; names[i++] = "r1";
      names[i++] = "s0"; names[i++] = "s1";
      names[i++] = "d0"; names[i++] = "d1";
      // not strictly part of the CG state, but used to initialize velocity solve
      names[i++] = "ux"; 
      if (SpaceDim > 1) names[i++] = "uy";		

      
      Vector<LevelData<FArrayBox>*> vdata(m_finest_level+1);
      for (int lev = 0; lev <= m_finest_level;lev++)
	{
	  vdata[lev] = new LevelData<FArrayBox>(m_grids[lev],names.size(),IntVect::Unit);
	  int j(0);
	  a_x[lev]->copyTo(Interval(0,1), *vdata[lev], Interval(j,j+1));j+=2;
	  a_r[lev]->copyTo(Interval(0,1), *vdata[lev], Interval(j,j+1));j+=2;
	  a_s[lev]->copyTo(Interval(0,1), *vdata[lev], Interval(j,j+1));j+=2;
	  a_d[lev]->copyTo(Interval(0,1), *vdata[lev], Interval(j,j+1));j+=2;
	  m_velb[lev]->copyTo(Interval(0,SpaceDim-1), *vdata[lev], Interval(j,j+SpaceDim-1));j+=SpaceDim;	
	}

      WriteAMRHierarchyHDF5(fileName, m_grids, vdata, names, m_domain[0].domainBox(),
			    m_dx[0][0], 0.0, m_time , m_refRatio, vdata.size());
      
    }
  
}

bool InverseVerticallyIntegratedVelocitySolver::readCGState
( int& a_iter,
  bool& a_fresh_start,
  Vector<LevelData<FArrayBox>* >& a_x,
  Vector<LevelData<FArrayBox>* >& a_r,
  Vector<LevelData<FArrayBox>* >& a_s,
  Vector<LevelData<FArrayBox>* >& a_d)
{

  if (m_config.m_CGreadStateFile == "")
    {
      //no file specified, so not reading...
      return false;
    }
  
  if (m_config.m_CGreadStateIter < a_iter)
    {
      //we read the file already...
      return false;
    }

  //
  pout() << " InverseVerticallyIntegratedVelocitySolver::readCGState: reading CG state from "
	 << m_config.m_CGreadStateFile << std::endl;

  a_fresh_start = m_config.m_CGfreshStart;
  
  // data we need to read an AMR hierarchy
  int numLevels;
  Vector<LevelData<FArrayBox>* > data;
  Vector<DisjointBoxLayout> grids;
  Vector<std::string> names;
  Vector<int> ratio;
  Real crseDx = 0.0, dt = 0.0, time = 0.0;
  Box crseBox;
  
  int status = ReadAMRHierarchyHDF5
    (m_config.m_CGreadStateFile,grids,data,names,crseBox,crseDx,dt,time,
     ratio,numLevels);

  CH_assert(data.size() == a_x.size()); 

  for (int lev = 0; lev < numLevels; lev++)
    {
      int j(0);
      data[lev]->copyTo(Interval(j,j+1), *a_x[lev], Interval(0,1));j+=2;
      data[lev]->copyTo(Interval(j,j+1), *a_r[lev], Interval(0,1));j+=2;
      data[lev]->copyTo(Interval(j,j+1), *a_s[lev], Interval(0,1));j+=2;
      data[lev]->copyTo(Interval(j,j+1), *a_d[lev], Interval(0,1));j+=2;
      if ((names.size() > j) && names[j] == "ux") 
	{
	   //save state contains the velocity
	   data[lev]->copyTo(Interval(j,j+SpaceDim-1), *m_velb[lev], Interval(0,SpaceDim-1));j+=SpaceDim;
	}
    }

  a_iter = m_config.m_CGreadStateIter;
  
  return true;
  
}



void InverseVerticallyIntegratedVelocitySolver::writeState
(const std::string& a_file, int a_counter,
 const Vector<LevelData<FArrayBox>* >& a_x,
 const Vector<LevelData<FArrayBox>* >& a_g) const
{

  
  
  pout() << "writing state to " << a_file << std::endl;
  Vector<std::string> names;
  names.resize(0);
  names.push_back("X0");
  names.push_back("X1");
  names.push_back("C");
  names.push_back("Cwshelf");
  names.push_back("muCoef");
  names.push_back("xVelb");
  names.push_back("yVelb");
  names.push_back("xVels");
  names.push_back("yVels");
  names.push_back("xVelo");
  names.push_back("yVelo");
  names.push_back("xAdjVel");
  names.push_back("yAdjVel");
  names.push_back("xAdjRhs");
  names.push_back("yAdjRhs");
  names.push_back("gradJC");
  names.push_back("gradJMuCoef");
  names.push_back("velc");
  names.push_back("thickness");
  names.push_back("Z_base");
  names.push_back("Z_surface");
  Vector<LevelData<FArrayBox>*> vdata(m_finest_level+1);
  for (int lev = 0; lev <= m_finest_level;lev++)
    {
      vdata[lev] = new LevelData<FArrayBox>(m_grids[lev],names.size(),IntVect::Zero);
      LevelData<FArrayBox>& data = *vdata[lev];
      int j = 0;
      a_x[lev]->copyTo(Interval(0,0),data,Interval(j,j));j++;
      a_x[lev]->copyTo(Interval(1,1),data,Interval(j,j));j++;
      m_C[lev]->copyTo(Interval(0,0),data,Interval(j,j));j++; 
      m_Cmasked[lev]->copyTo(Interval(0,0),data,Interval(j,j));j++; 
      m_muCoef[lev]->copyTo(Interval(0,0),data,Interval(j,j));j++; 
      m_velb[lev]->copyTo(Interval(0,1),data,Interval(j,j+1));j+=2;
      m_vels[lev]->copyTo(Interval(0,1),data,Interval(j,j+1));j+=2;
      m_velObs[lev]->copyTo(Interval(0,1),data,Interval(j,j+1));j+=2;
      m_adjVel[lev]->copyTo(Interval(0,1),data,Interval(j,j+1));j+=2;
      m_adjRhs[lev]->copyTo(Interval(0,1),data,Interval(j,j+1));j+=2;
      a_g[lev]->copyTo(Interval(0,0),data,Interval(j,j));j++;
      a_g[lev]->copyTo(Interval(1,1),data,Interval(j,j));j++;
      m_velCoef[lev]->copyTo(Interval(0,0),data,Interval(j,j));j++;
      m_coordSys[lev]->getH().copyTo(Interval(0,0),data,Interval(j,j));j++;
      m_coordSys[lev]->getTopography().copyTo(Interval(0,0),data,Interval(j,j));j++;
      m_coordSys[lev]->getSurfaceHeight().copyTo(Interval(0,0),data,Interval(j,j));j++;
    
    }
  const Real dt = 1.0;
  const Real time = Real(a_counter);
 
  for (int lev = vdata.size()-1;lev > 0 ; lev--)
    {
      LevelData<FArrayBox>& fine = *vdata[lev];
      LevelData<FArrayBox>& crse = *vdata[lev-1];
      CoarseAverage ca(m_grids[lev],fine.nComp(),m_refRatio[lev-1]);
      ca.averageToCoarse(crse,fine);
    }
 
  m_amrIce->writeAMRHierarchyHDF5(a_file,m_grids,vdata,names, m_domain[0].domainBox(),
   				  m_dx[0][0], dt, time , m_refRatio, vdata.size());
  
  for (int lev = 0; lev <= m_finest_level;lev++)
    {
      delete vdata[lev];
    }
}

std::string InverseVerticallyIntegratedVelocitySolver::CGstateFile(int a_iter) const
{
  std::stringstream ss;
  ss << m_config.m_CGsaveStateFileNameBase;
  
  ss.width(2);ss.fill('0');ss << m_finest_level;
  ss.width(0); ss << "lev.";
    
  ss.width(6);ss.fill('0');ss << int(m_time/m_config.m_dtTypical);
  //ss.width(0); ss << "t.";
    
  ss.width(6);ss.fill('0');ss << a_iter;
  ss.width(0);ss << ".2d.hdf5";
  return ss.str(); 
}    


std::string InverseVerticallyIntegratedVelocitySolver::outerStateFile() const
{
  std::stringstream ss;
  ss << m_config.m_outerStepFileNameBase;
  

  ss.width(2);ss.fill('0');ss << m_finest_level;
  ss.width(0); ss << "lev.";
    
  ss.width(6);ss.fill('0');ss << int(m_time/m_config.m_dtTypical);
  //ss.width(0); ss << "t.";
    
  ss.width(6);ss.fill('0');ss << m_outerCounter;
  ss.width(0);ss << ".2d.hdf5";
  return ss.str(); 
}    

std::string InverseVerticallyIntegratedVelocitySolver::innerStateFile() const
{
  std::stringstream ss;
  ss << m_config.m_innerStepFileNameBase;
  ss.width(2);ss.fill('0');ss << m_finest_level;
  ss.width(0); ss << "lev.";

  ss.width(6);ss.fill('0');ss << int(m_time*12.0);
  //ss.width(0); ss << "t.";
    
  ss.width(6);ss.fill('0');ss << m_innerCounter;
  ss.width(0);ss << ".2d.hdf5";



  return ss.str(); 
}

#include "NamespaceFooter.H"
