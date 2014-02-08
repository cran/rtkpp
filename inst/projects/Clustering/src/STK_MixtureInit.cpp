/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2013  Serge Iovleff

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : S..._DOT_I..._AT_stkpp.org (see copyright for ...)
*/

/*
 * Project:  stkpp::Clustering
 * created on: 24 août 2013
 * Author:   iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_MixtureInit.cpp
 *  @brief In this file we implement the initialization methods.
 **/

#include "STKernel/include/STK_Exceptions.h"
#include "../include/STK_MixtureInit.h"
#include "../include/STK_MixtureAlgo.h"
#include "../include/STK_IMixtureComposer.h"

namespace STK
{

/** destructor */
IMixtureInit::~IMixtureInit()
{ if (p_initAlgo_) delete p_initAlgo_;}

/* launch the algorithm of the initialization step */
bool IMixtureInit::runInitAlgo()
{
  if (p_initAlgo_)
  {
    p_initAlgo_->setModel(p_model_);
    return (p_initAlgo_->run());
  }
  msg_error_ = _T("p_initAlgo is not initialized.");
  return false;
}

/* call the randomClassInit() model initialization.
 * @return @c true if no error occur, @c false otherwise*/
bool RandomInit::run()
{
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("-------------------------------\n")
           << _T("Entering RandomInit::run()\n");
#endif
  bool result = false;
  for (int iTry= 0; iTry < nbTry_; ++iTry)
  {
     try
     {
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("p_model_->randomInit() starting.\n");
#endif
       p_model_->randomInit();
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("p_model_->randomInit() terminated.\n");
#endif
       if (runInitAlgo()) { result = true; break;}
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("runInitAlgo() failed\n");
#endif
     }
     catch (Clust::exceptions const& error)
     {
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("In RandomInit::run(), try number") << iTry << " generate an exception.\n";
#endif
     }
  } // iTry
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("Exiting RandomInit::run()\n")
           << _T("-------------------------------\n");
#endif
  if (result) return true;
  msg_error_ = _T("All initialization failed.");
  return false;
}

/* call the randomClassInit() model initialization.
 * @return @c true if no error occur, @c false otherwise*/
bool ClassInit::run()
{
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("-------------------------------\n")
           << _T("Entering ClassInit::run()\n");
#endif
  bool result = false;
  for (int iTry= 0; iTry < nbTry_; ++iTry)
  {
#ifdef STK_MIXTURE_VERBOSE
    stk_cout << _T("iTry = ") << iTry << "\n";
#endif
    try
    {
      p_model_->randomClassInit();
      if (runInitAlgo()) { result = true; break;}
    }
    catch (Clust::exceptions const& error)
    {
#ifdef STK_MIXTURE_VERBOSE
      stk_cout << _T("In ClassInit::run(), try number: ") << iTry << " generate exception.\n";
#endif
    }
  } // iTry
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("Exiting ClassInit::run()\n")
           << _T("-------------------------------\n");
#endif
  if (result) return true;
  // otherwise return false
  msg_error_ = _T("All initialization failed.");
  return false;
}

/* call the randomClassInit() model initialization.
 * @return @c true if no error occur, @c false otherwise*/
bool FuzzyInit::run()
{
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("-------------------------------\n")
           << _T("Entering FuzzyInit::run()\n");
#endif
  bool result = false;
  for (int iTry= 0; iTry < nbTry_; ++iTry)
  {
    try
    {
      p_model_->randomFuzzyInit();
      if (runInitAlgo()) { result = true; break;}
    }
    catch (Clust::exceptions const& error)
    {/* do nothing and retry*/
#ifdef STK_MIXTURE_VERY_VERBOSE
      stk_cout << _T("In FuzzyInit::run try:") << iTry << " failed.\n";
#endif
    }
  } // iTry
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("Exiting FuzzyInit::run()\n")
           << _T("-------------------------------\n");
#endif
  if (result) return true;
  // otherwise return false
  msg_error_ = _T("All initialization failed.");
  return false;
}

} // namespace STK
