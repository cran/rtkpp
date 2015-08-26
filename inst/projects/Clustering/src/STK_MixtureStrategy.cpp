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
 * Project:  stkpp::
 * created on: 3 sept. 2013
 * Author:   iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_MixtureStrategy.cpp
 *  @brief In this file we implement the strategies for estimating mixture model.
 **/

#include "STKernel/include/STK_Exceptions.h"
#include "../include/STK_MixtureStrategy.h"
#include "../include/STK_MixtureInit.h"
#include "../include/STK_MixtureAlgo.h"
#include "../include/STK_IMixtureComposer.h"

namespace STK
{

/* copy constructor
 *  @param strategy the strategy to copy
 **/
IMixtureStrategy::IMixtureStrategy( IMixtureStrategy const& strategy)
                                  : IRunnerBase(strategy), nbTry_(strategy.nbTry_)
                                  , p_model_(strategy.p_model_)
                                  , p_init_(strategy.p_init_->clone())
{}

/* store a model in p_model_ if it is better.
 * @param p_otherModel the model to store
 **/
void IMixtureStrategy::storeModel(IMixtureComposer*& p_otherModel)
{ if (p_model_->lnLikelihood()<p_otherModel->lnLikelihood())
  { std::swap(p_model_, p_otherModel);}
}

/* destructor */
IMixtureStrategy::~IMixtureStrategy() { if (p_init_) delete p_init_;}

/* destructor */
SimpleStrategyParam::~SimpleStrategyParam()
{ if (p_algo_) delete p_algo_;}

/* destructor */
XemStrategyParam::~XemStrategyParam()
{
  if (p_shortAlgo_) delete p_shortAlgo_;
  if (p_longAlgo_) delete p_longAlgo_;
}

/* destructor */
FullStrategyParam::~FullStrategyParam()
{
  if (p_shortAlgo_) delete p_shortAlgo_;
  if (p_longAlgo_) delete p_longAlgo_;
}

/* run the simple strategy */
bool SimpleStrategy::run()
{
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("-----------------------------------------------\n");
  stk_cout << _T("Entering SimpleStrategy::run() with: ")
           << _T("nbTry_ = ") << nbTry_ << _T("\n");
#endif
  IMixtureComposer* p_currentModel = 0;
  if (p_model_->state() < 1) { p_model_->randomFuzzyInit();}
  Real value = p_model_->lnLikelihood();
  try
  {
    p_currentModel = p_model_->create();

    // initialize and run algo. break if success
    for (int iTry = 0; iTry < nbTry_; ++iTry)
    {
      // initialize current model
      p_init_->setModel(p_currentModel);
      if (!p_init_->run())
      {
#ifdef STK_MIXTURE_VERBOSE
        stk_cout << "SimpleStrategy::run(), Initialization failed.\n";
#endif
        msg_error_ += STKERROR_NO_ARG(SimpleStrategy::run,Initialization failed\n);
        msg_error_ += p_init_->error();
      } // init step
      p_param_->p_algo_->setModel(p_currentModel);
      // initial
      if (!p_param_->p_algo_->run())
      {
#ifdef STK_MIXTURE_VERBOSE
        stk_cout << "SimpleStrategy::run(), Long run failed.\n";
#endif
        msg_error_ += STKERROR_NO_ARG(SimpleStrategy::run,long run failed\n);
        msg_error_ += p_param_->p_algo_->error();
      }
      else { break;}  // we get a result
    } // iTry
  }
  catch (Exception const& e)
  {
    msg_error_ = e.error();
    return false;
  }
  storeModel(p_currentModel);
  delete p_currentModel;
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << "SimpleStrategy::run() terminated.              \n";
  stk_cout << "-----------------------------------------------\n";
#endif
  if (p_model_->lnLikelihood() <= value)
  {
    msg_error_ += STKERROR_NO_ARG(SimpleStrategy::run,No gain\n);
    return false;
  }
  return true;
}

/* run the xem strategy */
bool XemStrategy::run()
{
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("Entering XemStrategy::run() with:\n")
           << _T("nbTry_ = ") << nbTry_ << _T("\n")
           << _T("nbShortRun_ = ") << p_param_->nbShortRun_ << _T("\n");
#endif
  // the current model is used in the short runs
  IMixtureComposer* p_currentModel     = p_model_->create();
  IMixtureComposer* p_currentBestModel = p_model_->create();
  // add some perturbation to the tik and compute the ln-likelihood
  if (p_model_->state() < 1) { p_model_->randomFuzzyInit();}
  Real value = p_model_->lnLikelihood();
  // start estimation
  try
  {
    for (int iTry = 0; iTry < nbTry_; ++iTry)
    {
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("-------------------------------\n")
           << _T("try number = ") << iTry << _T("\n");
#endif
      // find best of the shortModel and save it in p_currentBestModel
      for (int iShortRun = 0; iShortRun < p_param_->nbShortRun_; ++iShortRun)
      {
        // initialize current model
        p_init_->setModel(p_currentModel);
        if (p_init_->run())
        {
          // perform short run on the current model
          p_param_->p_shortAlgo_->setModel(p_currentModel);
          p_param_->p_shortAlgo_->run();
          // if we get a better result, swap it with currentBestModel
          if( p_currentBestModel->lnLikelihood()<p_currentModel->lnLikelihood())
          { std::swap(p_currentModel, p_currentBestModel);}
        } // initialization
      } // iShortRun
      // in case nbShortRun_==0
      // try to initialize bestCurrentModel, otherwise go to a next try
      if (p_param_->nbShortRun_ == 0)
      {
        // initialize current model
        p_init_->setModel(p_currentBestModel);
        if (!p_init_->run())
        { continue; }// model not initialized, we go to the next trial
      }
#ifdef STK_MIXTURE_VERBOSE
      stk_cout << _T("iTry =") << iTry
               << _T(". In XemStrategy::run(), short run terminated. best model:\n");
      p_currentBestModel->writeParameters(stk_cout);
      stk_cout << _T("\n\n");
#endif
      // start a long run with the better model
      p_param_->p_longAlgo_->setModel(p_currentBestModel);
      if (p_param_->p_longAlgo_->run()) { storeModel(p_currentBestModel); break;}
#ifdef STK_MIXTURE_VERBOSE
            stk_cout << "In FullStrategy::run(), Long run Failed.\n";
#endif
    } // end iTry
    delete p_currentBestModel;
    delete p_currentModel;
  } catch (Exception const& e)
  {
    msg_error_ = e.error();
    return false;
  }
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << "XemStrategy::run() terminated.\n";
  stk_cout << "-------------------------------\n";
#endif
  if (p_model_->lnLikelihood() <= value)
  {
    msg_error_ += STKERROR_NO_ARG(In XemStrategy::run,No gain\n);
    return false;
  }
  return true;
}

/* run the full strategy */
bool FullStrategy::run()
{
  // add some perturbation to the tik and compute the ln-likelihood
  p_model_->setState(Clust::modelInitialized_);
#ifdef STK_MIXTURE_DEBUG
  p_model_->writeParameters(stk_cout);
#endif
  p_model_->randomFuzzyInit();
#ifdef STK_MIXTURE_DEBUG
  p_model_->writeParameters(stk_cout);
#endif
  Real initialValue = p_model_->lnLikelihood();
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("<+++\n");
  stk_cout << _T("Entering FullStrategy::run() with nbTry_ = ") << nbTry_
           << _T(", nbShortRun_ = ") << p_param_->nbShortRun_
           << _T(", p_model_->lnLikelihood() = ") << initialValue
           << _T("\n");
#endif
  IMixtureComposer* p_bestModel = 0;
  IMixtureComposer* p_bestShortModel   = 0;
  // start estimation
  try
  {
    // Main loop. If the Full strategy success in estimating a model, the
    // iterations are stopped and the best model find is stored in p_model_
    for (int iTry = 0; iTry < nbTry_; ++iTry)
    {
      // in case nbShortRun_==0: initialize directly p_bestShortModel
      if (p_param_->nbShortRun_ <= 0)
      {
        if (!initStep(p_bestShortModel))
        {
          msg_error_ += STKERROR_NO_ARG(FullStrategy::run,init step failed\n);
          msg_error_ += p_param_->p_shortAlgo_->error();
#ifdef STK_MIXTURE_VERBOSE
          stk_cout << _T("In FullStrategy::run()") << _T(", iTyry =") << iTry
                   << _T(", init step failed\n");
          stk_cout << msg_error_ << _T("\n");
#endif
        }
      }
      else
      {
#ifdef STK_MIXTURE_VERY_VERBOSE
        stk_cout << _T("In FullStrategy::run(), entering short run steps\n")
                 << _T("iTyry =") << iTry << _T("\n");
#endif
        Real valueBest = -Arithmetic<Real>::infinity();
        for (int iShort=0; iShort < p_param_->nbShortRun_; ++iShort)
        {
          // perform nbInitRun_ initialization step and get the best result in p_bestModel
          if (!initStep(p_bestModel))
          {
            msg_error_ += STKERROR_NO_ARG(FullStrategy::run,init step failed\n);
            msg_error_ += p_param_->p_shortAlgo_->error();
#ifdef STK_MIXTURE_VERBOSE
            stk_cout << _T("In FullStrategy::run()") << _T(", iTyry =") << iTry << _T(", iShort =") << iShort
                     << _T(", init step failed\n");
            stk_cout << msg_error_ << _T("\n");
#endif
          }
          // In case an error occur in initStep
          if (!p_bestModel)
          {
            p_bestModel = p_model_->clone();
          }
          // perform short run with the current best model
          p_param_->p_shortAlgo_->setModel(p_bestModel);
          if (!p_param_->p_shortAlgo_->run())
          {
            msg_error_ += STKERROR_NO_ARG(FullStrategy::run,short algo failed\n);
            msg_error_ += p_param_->p_shortAlgo_->error();
#ifdef STK_MIXTURE_VERBOSE
            stk_cout << _T("In FullStrategy::run()") << _T(", iTyry =") << iTry << _T(", iShort =") << iShort
                     << _T(", short Algo fail\n");
            stk_cout << msg_error_ << _T("\n");
#endif
          }
          // if we get a better result, store it in p_bestShortModel
          Real value = p_bestModel->lnLikelihood();
          if( valueBest<value)
          {
            std::swap(p_bestShortModel, p_bestModel);
            valueBest  = value;
#ifdef STK_MIXTURE_VERY_VERBOSE
            stk_cout << _T("In FullStrategy::run()")
                     << _T(", iTyry =") << iTry << _T(", iShort =") << iShort
                     << _T(", get better value in short run. valueBest =") << valueBest << _T("\n");
#endif
          }
        } // ishort
        // release memory
        if (p_bestModel)
        {
          delete p_bestModel; p_bestModel = 0;
        }
      }
      // in case all initialization failed
      if (!p_bestShortModel){ p_bestShortModel = p_model_->clone();}
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("In FullStrategy::run() all short run") << _T(", iTyry =") << iTry  << _T(" terminated.\n")
           << _T("p_bestShortModel->lnLikelihood() = ") << p_bestShortModel->lnLikelihood()
           << _T("\n");
#endif
      // start a long run with p_bestShortModel. If success, save model
      // and exit the iTry loop
      p_param_->p_longAlgo_->setModel(p_bestShortModel);
      if (!p_param_->p_longAlgo_->run())
      {
        msg_error_ += STKERROR_NO_ARG(FullStrategy::run,long algo failed\n);
        msg_error_ += p_param_->p_longAlgo_->error();
#ifdef STK_MIXTURE_VERBOSE
        stk_cout << _T("In FullStrategy::run(): Long Algo failed\n");
#endif
      }
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("In FullStrategy::run() long run") << _T(", iTyry =") << iTry  << _T(" terminated.\n")
           << _T("p_bestShortModel->lnLikelihood() = ") << p_bestShortModel->lnLikelihood()
           << _T("\n");
#endif
      // if we get a better result, store it in p_model_ and stop to try
      if( p_model_->lnLikelihood()<p_bestShortModel->lnLikelihood())
      {
        std::swap(p_model_, p_bestShortModel);
        break;
      }
#ifdef STK_MIXTURE_VERBOSE
      stk_cout << _T("In FullStrategy::run(), iTry =") << iTry << _T(" failed\n");
#endif
      // release memory before next try
      if (p_bestModel)      delete p_bestModel; p_bestModel = 0;
      if (p_bestShortModel) delete p_bestShortModel; p_bestShortModel = 0;
    } // end iTry
  }
  catch (Exception const& e)
  {
    if (p_bestModel)      delete p_bestModel;
    if (p_bestShortModel) delete p_bestShortModel;
    msg_error_ += e.error();
    return false;
  }
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << "FullStrategy::run() terminated. \n";
  stk_cout << _T("+++>\n");
#endif
  // normally not needed
  if (p_bestModel) delete p_bestModel;
  if (p_bestShortModel)   delete p_bestShortModel;
  if (p_model_->lnLikelihood() <= initialValue)
  {
    msg_error_ += STKERROR_NO_ARG(FullStrategy::run,No gain\n);
    return false;
  }
  return true;
}

/* Perform the Initialization step*/
bool FullStrategy::initStep(IMixtureComposer*& p_bestModel)
{
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("<+\n");
  stk_cout << _T("Entering FullStrategy::initStep\n");
  stk_cout << _T("nbInitRun = ") <<  p_param_->nbInitRun_ << _T("\n");
#endif
  IMixtureComposer* p_initModel = p_model_->create();
  try
  {
    Real valueBest = -Arithmetic<Real>::infinity();
    for (int iInitRun=0; iInitRun < p_param_->nbInitRun_; iInitRun++)
    {
      // set current model
      p_init_->setModel(p_initModel);
      if (!p_init_->run())
      {
#ifdef STK_MIXTURE_VERBOSE
        stk_cout<< _T("FullStrategy::initStep, run failed:\n");
        stk_cout<< p_init_->error() << _T("\n");
#endif
      }
      Real value = p_initModel->lnLikelihood();
      // if we get a better result, swap it with currentBestModel
      if( valueBest<value)
      {
        std::swap(p_initModel, p_bestModel);
        valueBest = value;
#ifdef STK_MIXTURE_VERY_VERBOSE
      stk_cout << _T("FullStrategy::initStep, iInitRun ") << iInitRun
               << _T(", currentBest =") << valueBest << _T("\n");
#endif
        // in case p_bestModel was 0 pointer and there is more iterations
        if (!p_initModel && iInitRun <= p_param_->nbInitRun_)
        {
          p_initModel = p_model_->create();
        }
      }
    }
  }
  catch (Exception const& e)
  {
    // in case all initialization failed
    if (!p_bestModel)
    {
      p_bestModel = p_model_->clone();
    }
    if (p_initModel)
    {
      delete p_initModel; p_initModel = 0;
    }
    msg_error_ += e.error();
    return false;
  }
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("FullStrategy::initStep done\n");
  stk_cout << _T("p_bestModel->lnLikelihood() = ") <<  p_bestModel->lnLikelihood() << _T("\n");
  stk_cout << _T("+>\n");
#endif
  // in case all initialization failed or nbInitRun_ <= 0
  delete p_initModel; p_initModel = 0;

  if (!p_bestModel) p_bestModel = p_model_->clone();
  return true;
}

} // namespace STK



