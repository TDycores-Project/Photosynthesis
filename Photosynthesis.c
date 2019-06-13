static char help[] = "\nThe following implements photosynthesis as outlined in:\n\n\
A.P. Walker et al: Multi-assumption architecture and testbed (MAAT v.10),\n\
     Geosci. Model Dev., 11, 3159-3185, 2018.\n\n";
#include "tdycore.h"
#include "time.h"

typedef struct {
  PetscReal Rd,Vcmax,alphaT,Kc,Ko,O,ko,kc,TPU,tau,Gamma;
  PetscReal a,alphai,II,Jmax,J;
  PetscReal p,Ca,ri,g0,g1m,D,tb,U,d1,kappar,rb;
} TDyPhotoParams;

PetscReal _r(){
  return ((PetscReal)rand())/((PetscReal)RAND_MAX);
}

PetscErrorCode TDyPhotoDefaultParams(PetscBag bag){
  PetscErrorCode ierr;
  PetscFunctionBegin;

  TDyPhotoParams *params;
  ierr = PetscBagGetData(bag,(void **)&params);
  
  // A1: Carbon assimilation
  
  ierr = PetscBagRegisterReal(bag,&params->Rd    ,_r()        ,"Rd","non-photorespiratory respiration CO2 [1e-6 mol m-2 s-1]"); CHKERRQ(ierr);
  ierr = PetscBagRegisterReal(bag,&params->Vcmax ,_r()        ,"Vcmax","CO2 [1e-6 mol m-2 s-1]"); CHKERRQ(ierr);
  ierr = PetscBagRegisterReal(bag,&params->alphaT,_r()        ,"alphaT","fraction of exported triose phosphate not returned [1]"); CHKERRQ(ierr);
  ierr = PetscBagRegisterReal(bag,&params->Kc    ,_r()        ,"Kc","Michaelis-Menten constants of RuBisCO for CO2 [Pa]"); CHKERRQ(ierr);
  ierr = PetscBagRegisterReal(bag,&params->Ko    ,_r()        ,"Ko","Michaelis-Menten constants of RuBisCO for O2 [1e3 Pa]"); CHKERRQ(ierr);
  ierr = PetscBagRegisterReal(bag,&params->O     ,_r()        ,"O","chloroplast O2 partial pressure, assumed atmospheric [1e3 Pa]"); CHKERRQ(ierr);
  ierr = PetscBagRegisterReal(bag,&params->ko    ,_r()        ,"ko","turnover rate of RuBisCO for oxygenation [s-1]"); CHKERRQ(ierr);
  ierr = PetscBagRegisterReal(bag,&params->kc    ,_r()        ,"kc","turnover rate of RuBisCO for carboxylation [s-1]"); CHKERRQ(ierr);
  ierr = PetscBagRegisterReal(bag,&params->TPU   ,_r()        ,"TPU","Triose phosphate utilization CO2 [1e-6 mol m-2 s-1]"); CHKERRQ(ierr);
  ierr = PetscBagRegisterReal(bag,&params->tau   ,_r()        ,"tau"," *** CO2-O2 specificity ratio of RuBisCO"); CHKERRQ(ierr);
  ierr = PetscBagRegisterReal(bag,&params->Gamma ,_r()        ,"Gamma","*** photorespitory compensation point"); CHKERRQ(ierr);
  
  // A2: Electron transport

  ierr = PetscBagRegisterReal(bag,&params->a     ,_r()        ,"a","leaf absorptance"); CHKERRQ(ierr);
  ierr = PetscBagRegisterReal(bag,&params->alphai,_r()        ,"alphai","intrinsic quantum efficiency"); CHKERRQ(ierr);
  ierr = PetscBagRegisterReal(bag,&params->II    ,_r()        ,"II","incident photosynthetically active radiation [1e-6 mol m-2 s-1]"); CHKERRQ(ierr);
  ierr = PetscBagRegisterReal(bag,&params->Jmax  ,_r()        ,"Jmax","maximum electron transport rate"); CHKERRQ(ierr);
  ierr = PetscBagRegisterReal(bag,&params->J     ,_r()        ,"J","*** electron transport rate [ 1e-6 mol m-2 s-1]"); CHKERRQ(ierr);

  // A3: CO2 diffusion and resistance and stomatal conductance
  
  ierr = PetscBagRegisterReal(bag,&params->p     ,1.013250e-01,"p","atmospheric pressure [1e6 Pa]"); CHKERRQ(ierr);
  ierr = PetscBagRegisterReal(bag,&params->Ca    ,_r()        ,"Ca","atmospheric CO2 partial pressure [Pa]"); CHKERRQ(ierr);
  ierr = PetscBagRegisterReal(bag,&params->ri    ,_r()        ,"ri","internal resistance, CO2 [m2 s mol-1]"); CHKERRQ(ierr);
  ierr = PetscBagRegisterReal(bag,&params->g0    ,_r()        ,"g0","minimum stomatal conductance"); CHKERRQ(ierr);
  ierr = PetscBagRegisterReal(bag,&params->g1m   ,_r()        ,"g1m","stomatal conductance slope from Medlyn"); CHKERRQ(ierr);
  ierr = PetscBagRegisterReal(bag,&params->D     ,_r()        ,"D","vapor pressure deficit [1e3 Pa]"); CHKERRQ(ierr);
  ierr = PetscBagRegisterReal(bag,&params->tb    ,_r()        ,"tb","turbulent transfer coefficient [m s-0.5]"); CHKERRQ(ierr);
  ierr = PetscBagRegisterReal(bag,&params->U     ,_r()        ,"U","wind speed [m s-1]"); CHKERRQ(ierr);
  ierr = PetscBagRegisterReal(bag,&params->d1    ,_r()        ,"d1","leaf dimension in wind direction [m]"); CHKERRQ(ierr);
  ierr = PetscBagRegisterReal(bag,&params->kappar,_r()        ,"kappar","converts [s m-1] to [m2 s mol-1]"); CHKERRQ(ierr);
  ierr = PetscBagRegisterReal(bag,&params->rb    ,_r()        ,"rb","*** boundary layer resistance, H2O [m2 s mol-1]"); CHKERRQ(ierr);
    
  PetscFunctionReturn(0);
}

PetscErrorCode TDyUpdatePhotoParams(TDyPhotoParams *p)
{
  PetscFunctionBegin;
  p->rb = 1/p->tb/PetscSqrtReal(p->U/p->d1)*p->kappar; // boundary layer resistance, H2O [m2 s mol-1]
  p->J = p->a * p->alphai * p->II / PetscSqrtReal( 1 + PetscSqr( p->a * p->alphai * p->II / p->Jmax ) ); // electron transport rate [ 1e-6 mol m-2 s-1]
  p->tau = p->Ko * p->kc / ( p->Kc * p->ko ); // CO2-O2 specificity ratio of RuBisCO
  p->Gamma = 0.5 * p->O/ p->tau; // photorespitory compensation point
  PetscFunctionReturn(0);
}

/*
  Monolithic residual in terms of the net CO2 assimilation [1e-6 mol m-2 s-1]

  R[A] = 0

 */
PetscErrorCode ScalarResidual(SNES snes,Vec U,Vec R,void *ctx)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscScalar *rarray,A,Cb,Cbm,gs,rs,r,Cc,Acg,Ajg,Apg,Ag;
  const PetscScalar *uarray;
  TDyPhotoParams *p = (TDyPhotoParams *)ctx;
  ierr = VecGetArrayRead(U,&uarray);CHKERRQ(ierr);
  ierr = VecGetArray(R,&rarray);CHKERRQ(ierr);

  A = uarray[0];
  
  Cb = p->Ca - 1.4 * p->rb * A * p->p; // (A11a)
  Cbm = Cb / p->p; // Cb in molar units
  gs  = p->g0 + ( 1 + p->g1m / PetscSqrtReal(p->D) ) * A / Cbm; // stomatal conductance [mol m-2 s-1]  
  rs = 1/gs; // stomatal resistance, H2O [m2 s mol-1]
  r = 1.4 * p->rb + 1.6 * rs + p->ri; // resistance to CO2 diffusion [m2 s mol-1]
  Cc = p->Ca - r * A * p->p; // chloroplast CO2 partial pressure [Pa] (A9)
  Acg = p->Vcmax * Cc / (Cc + p->Kc * (1 + p->O / p->Ko) ); 
  Ajg = 0.25 * p->J * ( Cc / ( Cc + 2*p->Gamma ) );
  Apg = 3 * p->TPU * Cc / ( Cc + ( 1 + 3*p->alphaT ) * p->Gamma );
  Ag  = PetscMin(PetscMin(Acg,Ajg),Apg); // gross carboxylation rate CO2 [1e-6 mol m-2 s-1]

  rarray[0] = A - Ag * ( 1 - p->Gamma / Cc ) + p->Rd;
  
  ierr = VecRestoreArrayRead(U,&uarray);CHKERRQ(ierr);
  ierr = VecRestoreArray(R,&rarray);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


int main(int argc, char **argv)
{
  PetscErrorCode ierr;
  ierr = PetscInitialize(&argc,&argv,NULL,help);CHKERRQ(ierr);
  srand(time(0));
  
  /* Create parameters which are automatically given commandline options */
  PetscBag bag;
  ierr = PetscBagCreate(PETSC_COMM_WORLD,sizeof(TDyPhotoParams),&bag); CHKERRQ(ierr);
  ierr = PetscBagSetName(bag,"Default Photosynthesis Parameters:","those marked *** are computed from others"); CHKERRQ(ierr);
  ierr = TDyPhotoDefaultParams(bag); CHKERRQ(ierr);
  ierr = PetscBagSetFromOptions(bag); CHKERRQ(ierr);
  
  /* Update parameters before solve */
  TDyPhotoParams *params;
  ierr = PetscBagGetData(bag,(void **)&params);
  ierr = TDyUpdatePhotoParams(params); CHKERRQ(ierr);

  /* Solve the nonlinear system */
  SNES snes;
  ierr = SNESCreate(PETSC_COMM_SELF,&snes);CHKERRQ(ierr);
  ierr = SNESSetFunction(snes,NULL,ScalarResidual,params);CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  Vec U;
  ierr = VecCreate(PETSC_COMM_SELF,&U);CHKERRQ(ierr);
  ierr = VecSetType(U,VECSEQ);CHKERRQ(ierr);
  ierr = VecSetSizes(U,1,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = PetscBagView(bag,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ierr = SNESSolve(snes,NULL,U);CHKERRQ(ierr);
  
  ierr = VecDestroy(&U);CHKERRQ(ierr);
  ierr = SNESDestroy(&snes);CHKERRQ(ierr);  
  ierr = PetscBagDestroy(&bag); CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);
  return(0);
}

#ifdef __SAVE_FOR_LATER__
    /* Write parameters to a ASCII file (cannot read) */
  PetscViewer view;
  ierr = PetscViewerCreate(PETSC_COMM_WORLD,&view); CHKERRQ(ierr);
  ierr = PetscViewerSetType(view,PETSCVIEWERASCII); CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(view,FILE_MODE_WRITE); CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(view,"default.txt"); CHKERRQ(ierr);
  ierr = PetscBagView(bag,view); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&view); CHKERRQ(ierr);

  /* Write parameters to a binary file */
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"default.prm",FILE_MODE_WRITE,&view);
  ierr = PetscBagView(bag,view); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&view); CHKERRQ(ierr);

  /* Read parameters from a binary file */
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"default.prm",FILE_MODE_READ,&view);
  ierr = PetscBagLoad(view,bag); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&view); CHKERRQ(ierr);
#endif



