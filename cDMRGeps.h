//
// Based on dmrg.h  [https://itensor.org]
// Created by Erich Mueller in 2019
// Last edited by Shovan Dutta on 7/28/21
//

// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

#ifndef __ITENSOR_cDMRGeps_H
#define __ITENSOR_cDMRGeps_H

#include "itensor/iterativesolvers.h"
#include "itensor/mps/localmposet.h"
#include "itensor/mps/localmpo_mps.h"
#include "itensor/mps/sweeps.h"
#include "itensor/mps/DMRGObserver.h"
#include "itensor/util/cputime.h"
#include "localmpoproj.h"    // Functions of local Hamiltonian and discontinuity
#include "cio.h"     // Input-output functions (for saving observers on the fly)


namespace itensor {

    
template<class LocalOpT>
Real
cDMRGWorker(MPS & psi,                     // State
           LocalOpT & PH,                  // Hamiltonian and discontinuity
           Sweeps const& sweeps,           // Sweep parameters
           Args& args = Args::global());   // Other arguments

template<class LocalOpT>
Real
cDMRGWorker(MPS & psi,
           LocalOpT & PH,
           Sweeps const& sweeps,
           DMRGObserver & obs,             // Custom observers
           Args& args = Args::global());

//
// DMRG with given MPOs
//
Real inline
cDMRGeps(MPS & psi,
     MPO const& H,
     MPO const& penalty,
     Sweeps const& sweeps,
     Args& args = Args::global())
    {
    LocalMPOproj PH(H,penalty,args); //class of the two MPOs with defined functions
    Real energy = cDMRGWorker(psi,PH,sweeps,args); //active block
    return energy;
    }


//
// DMRGWorker
//
template<class LocalOpT>
Real
cDMRGWorker(MPS & psi,
           LocalOpT & PH,
           Sweeps const& sweeps,
           Args& args)
    {
    DMRGObserver obs(psi,args);
    Real energy = cDMRGWorker(psi,PH,sweeps,obs,args); //active block
    return energy;
    }


template<class LocalOpT>
Real
cDMRGWorker(MPS & psi,
           LocalOpT & PH,
           Sweeps const& sweeps,
           DMRGObserver & obs,
           Args& args)
    {
        println("In worker");
        
    if( args.defined("WriteM") )
      {
      if( args.defined("WriteDim") )
        {
        Global::warnDeprecated("Args WirteM and WriteDim are both defined. WriteM is deprecated in favor of WriteDim, WriteDim will be used.");
        }
      else
        {
        Global::warnDeprecated("Arg WriteM is deprecated in favor of WriteDim.");
        args.add("WriteDim",args.getInt("WriteM"));
        }
      }

    // Truncate blocks of degenerate singular values (or not)
    args.add("RespectDegenerate",args.getBool("RespectDegenerate",true));

    const bool silent = args.getBool("Silent",false);
    if(silent)
        {
        args.add("Quiet",true);
        args.add("PrintEigs",false);
        args.add("NoMeasure",true);
        args.add("DebugLevel",0);
        }
    const bool quiet = args.getBool("Quiet",false);
    const int debug_level = args.getInt("DebugLevel",(quiet ? 0 : 1));

    const int N = length(psi);
    Real Ewpenalty = NAN; //energy with penalty - updated after each step

    Real oldenergywpenalty = NAN; //energy with penalty - updated after each sweep
    Real soldenergywpenalty = NAN; //temporary variable for energy update
        
    psi.position(1); //right normalize

    args.add("DebugLevel",debug_level);
    args.add("DoNormalize",true);
        
        //
        // Load new arguments
        //
        auto MaxSweeps=args.getInt("MaxSweeps",sweeps.nsweep()); //max sweep number
        if (MaxSweeps>sweeps.nsweep())
        {
            MaxSweeps=sweeps.nsweep();
        }
        auto MinSweeps = args.getInt("MinSweeps",1); //min sweep number
        auto SweepThresh = args.getReal("SweepThresh",1E-8); //convergence threshold
        auto eps=args.getReal("eps",0.01); //eps = 1/penalty
        auto savethresh=args.getReal("SaveThresh",1E-6); //truncate saved numbers
        auto storelocaleverysweep=args.getBool("StoreLocalHistory",false); //dms
        auto localsweepname=args.getString("LocalSweepName"); //file name for dms
        auto saveenergy=args.getBool("saveenergyhistory",false); //energy each step
        auto energyname=args.getString("ehistname"); //file name for energy
    
        //
        // Open files for outputting intermediate results
        //
        std::ofstream efile; //energy file
        std::ofstream lsweepfile; //local dms file
        
        println("Open output files");
        
        if (saveenergy)
        {
            efile.open(energyname,std::ofstream::app); // open to append
            efile.precision(16); //16-digit precision
            efile.setf(std::ios::fixed);
        }

        if (storelocaleverysweep)
        {
            lsweepfile.open(localsweepname,std::ofstream::app); // open to append
            lsweepfile.precision(16);
            lsweepfile.setf(std::ios::fixed);
        }
        
        
        // New control structure
        bool notconverged=true;
        
        std::cout<<"starting loop, max="<<sweeps.nsweep()<<"\n";
        
        // Set weight of the discontinuity operator
        PH.weight(1/eps);
    
    
    for(int sw = 1; (sw <= sweeps.nsweep()) and notconverged; ++sw) //start sweeps
        {
        std::cout<<"Starting Sweep number "<<sw;
            
        cpu_time sw_time; //measure sweep durations //alternatively use clock()
        args.add("Sweep",sw);
        args.add("NSweep",sweeps.nsweep());
        args.add("Cutoff",sweeps.cutoff(sw));
        args.add("MinDim",sweeps.mindim(sw));
        args.add("MaxDim",sweeps.maxdim(sw));
        args.add("Noise",sweeps.noise(sw));
        args.add("MaxIter",sweeps.niter(sw));



        for(int b = 1, ha = 1; ha <= 2; sweepnext(b,ha,N)) //sweep over sites b
            {
            if(!quiet)
                {
                printfln("Sweep=%d, HS=%d, Bond=%d/%d",sw,ha,b,(N-1));
                }

            PH.position(b,psi); //orthogonality center at b

            auto phi = psi(b)*psi(b+1); //two-site optimization
                
                std::cout <<"Calling Davidson with niter="<<args.getInt("MaxIter")<<"\n";
                Ewpenalty=davidson(PH,phi,args); //minimization step
                auto energy=PH.Hexpect(phi); //final energy
                auto discontinuity=PH.Projexpect(phi); //final discontinuity
            
            auto spec = psi.svdBond(b,phi,(ha==1?Fromleft:Fromright),PH,args); //split bond by SVD
                
            if(!quiet)
                {
                    std::cout<<"\n";
                    std::cout<<"Ran Davidson with discontinuity weights: "<<1/eps<<"\n";
                    std::cout<<"Davidson Outputs: "<<Ewpenalty<<"\n";
                    std::cout<<"Energy: "<<energy<<"\n";
                    std::cout<<"Discontinuity: "<<discontinuity<<"\n";
                    std::cout<<"\n";
                printfln("    Truncated to Cutoff=%.1E, Min_dim=%d, Max_dim=%d",
                          sweeps.cutoff(sw),
                          sweeps.mindim(sw), 
                          sweeps.maxdim(sw) );
                printfln("    Trunc. err=%.1E, States kept: %s",
                         spec.truncerr(),
                         showDim(linkIndex(psi,b)) );
                }
                
               

                
            obs.lastSpectrum(spec);

            args.add("AtBond",b);
            args.add("HalfSweep",ha);
            args.add("Ewpenalty",Ewpenalty);
            args.add("energy",energy);
            args.add("discontinuity",discontinuity);
            args.add("Truncerr",spec.truncerr()); 

            obs.measure(args);
                
                auto sm = sw_time.sincemark(); //CPU and wall times for this step
                
                // Save intermediate results
                if (saveenergy)
                {
                    efile << ",{\"Ewpenalty\"->"<<Ewpenalty<<",\"disc\"->"<<discontinuity;
                    efile << ",\"E\"->"<<energy;
                    efile<< ",\"center\"->"<<b<<",\"dir\"->"<<ha;
                    efile<<",\"cputime\"->"<<sm.time<<",\"walltime\"->"<<sm.wall;
                    efile<< ",\"bonddim\"->"<<dim(linkIndex(psi,b))<<",\"eps\"->"<<eps<<"}";
                    
                }
                

            } //for loop over sites b

            
        if(!silent)
            {
            auto sm = sw_time.sincemark(); //sweep duration
            printfln("    Sweep %d/%d CPU time = %s (Wall time = %s)",
                      sw,sweeps.nsweep(),showtime(sm.time),showtime(sm.wall));
            }

        
            // Save density mats after each sweep
            if (storelocaleverysweep)
            {
                lsweepfile<<",";
                savelocal(lsweepfile,psi,savethresh); //saveimg
            }
            
            
            if(obs.checkDone(args)) break;
    
            if((sw>MinSweeps) and (abs(oldenergywpenalty/Ewpenalty-1)<SweepThresh)) //check for convergence
            {
                notconverged=false;
                std::cout<<"\nAccuracy threshold satisfied after "<<sw<<" sweeps\n";
                std::cout<<"MinSweeps= "<<MinSweeps<<" MaxSweeps="<<MaxSweeps<<" Threshold="<<SweepThresh<<"\n";
                std::cout<<"Last two energies with penalties: "<<Ewpenalty<<","<<oldenergywpenalty<<"\n\n";
                
            }
            // Update old and new energies
            soldenergywpenalty=oldenergywpenalty;
            oldenergywpenalty=Ewpenalty;
            
            args.add("numsweeps",sw);
            
            
        } //for loop over sweeps
        
    
            if (notconverged) //print if convergence threshold was not met
            {
                std::cout<<"\nAccuracy threshold "<<SweepThresh<<" not satisfied after "<<MaxSweeps<<" sweeps.\n";
                std::cout<<"Last two energies with penalties: "<<Ewpenalty<<","<<soldenergywpenalty<<" Difference:"<<soldenergywpenalty/Ewpenalty-1<<"\n\n";
            }

        
        // Close Files
        if (saveenergy)
        {
            efile.close();
        }
        
        if (storelocaleverysweep)
        {
            lsweepfile.close();
        }
        
        
    psi.normalize();

    return Ewpenalty; //final energy with penalty
    }


} //namespace itensor


#endif
