/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * timing.c
 *
 * This file contains routines that deal with timing Metis
 *
 * Started 7/24/97
 * George
 *
 * $Id: timing.c,v 1.2 2002/08/10 06:29:34 karypis Exp $
 *
 */

#include <metislib.h>


/*************************************************************************
* This function clears the timers
**************************************************************************/
void InitTimers(CtrlType *ctrl)
{
  gk_clearcputimer(ctrl->TotalTmr);
  gk_clearcputimer(ctrl->InitPartTmr);
  gk_clearcputimer(ctrl->MatchTmr);
  gk_clearcputimer(ctrl->ContractTmr);
  gk_clearcputimer(ctrl->CoarsenTmr);
  gk_clearcputimer(ctrl->UncoarsenTmr);
  gk_clearcputimer(ctrl->RefTmr);
  gk_clearcputimer(ctrl->ProjectTmr);
  gk_clearcputimer(ctrl->SplitTmr);
  gk_clearcputimer(ctrl->SepTmr);
  gk_clearcputimer(ctrl->AuxTmr1);
  gk_clearcputimer(ctrl->AuxTmr2);
  gk_clearcputimer(ctrl->AuxTmr3);
  gk_clearcputimer(ctrl->AuxTmr4);
  gk_clearcputimer(ctrl->AuxTmr5);
  gk_clearcputimer(ctrl->AuxTmr6);
}



/*************************************************************************
* This function prints the various timers
**************************************************************************/
void PrintTimers(CtrlType *ctrl)
{
  mprintf("\nTiming Information -------------------------------------------------");
  mprintf("\n Multilevel: \t\t %7.3f", gk_getcputimer(ctrl->TotalTmr));
  mprintf("\n     Coarsening: \t\t %7.3f", gk_getcputimer(ctrl->CoarsenTmr));
  mprintf("\n            Matching: \t\t\t %7.3f", gk_getcputimer(ctrl->MatchTmr));
  mprintf("\n            Contract: \t\t\t %7.3f", gk_getcputimer(ctrl->ContractTmr));
  mprintf("\n     Initial Partition: \t %7.3f", gk_getcputimer(ctrl->InitPartTmr));
  mprintf("\n   Construct Separator: \t %7.3f", gk_getcputimer(ctrl->SepTmr));
  mprintf("\n     Uncoarsening: \t\t %7.3f", gk_getcputimer(ctrl->UncoarsenTmr));
  mprintf("\n          Refinement: \t\t\t %7.3f", gk_getcputimer(ctrl->RefTmr));
  mprintf("\n          Projection: \t\t\t %7.3f", gk_getcputimer(ctrl->ProjectTmr));
  mprintf("\n     Splitting: \t\t %7.3f", gk_getcputimer(ctrl->SplitTmr));
  mprintf("\n          AUX1: \t\t %7.3f", gk_getcputimer(ctrl->AuxTmr1));
  mprintf("\n          AUX2: \t\t %7.3f", gk_getcputimer(ctrl->AuxTmr2));
  mprintf("\n          AUX3: \t\t %7.3f", gk_getcputimer(ctrl->AuxTmr3));
  mprintf("\n********************************************************************\n");
}



