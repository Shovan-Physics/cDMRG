//
//  cdmrg.cpp
//  
//  ver 1.01
//  Created by Erich Mueller on 6/16/19.
//  Edited 7/30/19
//  Edited 3/2/20
//  Edited by Shovan Dutta on 5/6/20
//  Edited 7/28/21


#include <stdio.h>
#include "itensor/all.h"   // Note: cpu_time defined in "itensor/util/cputime.h"
#include "cSite.h"         // Site set
#include "cDMRGeps.h"      // DMRG sweeps for a given discontinuity penalty (1/eps)
#include "cio.h"           // Input-output functions
#include "readvec.h"       // Read vector inputs
#include <fstream>
#include "localmpoproj.h"  // Operations with local Hamiltonian and discontinuity


using namespace std;
using namespace itensor;


inline bool fileexists (const std::string& name) // Check if a file exists
{
    return (access(name.c_str(), F_OK) != -1);
}


int main(int argc, char* argv[])  // Main module
{
    
    if(argc < 2)
    {
        printfln("Usage: %s input_file",argv[0]);
        println("");
        println("Input file can be created using Mathematica interface");
        return 0;
    }
    
    cpu_time overall_time; // Measure overall CPU- and wall times
    
    // Alternative way to measure times - used in benchmarking:
    // std::clock_t c_i = std::clock(); //initial CPU time
    // auto t_i = std::chrono::high_resolution_clock::now(); //initial wall time
    
    
//----------------------Load parameters from input file---------------------------//
    
    
    std::string inputfile(argv[1]);
    auto rawname = inputfile.substr(0, inputfile.find_last_of("."));
    println("Loading parameters from ",rawname);
    
    auto parameters = InputGroup(argv[1],"parameters"); //get parameters
    auto numsites = parameters.getInt("numsites"); //number of segments
    auto loadfromfile = parameters.getYesNo("loadfromfile",false); //initial state
    auto targetdiscontinuity = parameters.getReal("targetdiscontinuity",1.e-12);
    auto savethresh = parameters.getReal("savethresh",0.0); //truncate saved numbers
    
    // Parameters for DMRG cycles with different penalty (1/eps)
    auto numeps = parameters.getInt("numeps",0); //how many cycles
    auto epsgroup = InputGroup(argv[1],"eps");
    auto epsvals = ReadRealVec(epsgroup,numeps); //eps values (=1/penalty)
    auto MaxItergroup = InputGroup(argv[1],"MaxIters");
    auto MaxItervals = ReadIntVec(MaxItergroup,numeps); //max eigensolver iterations
    auto cutoffgroup = InputGroup(argv[1],"cutoffs");
    auto cutoffvals = ReadRealVec(cutoffgroup,numeps); //singular-value cutoffs
    auto MaxDimgroup = InputGroup(argv[1],"MaxDim");
    auto MaxDimvals = ReadIntVec(MaxDimgroup,numeps); //max bond dimensions
    auto threshgroup = InputGroup(argv[1],"epsthresh");
    auto threshvals = ReadRealVec(threshgroup,numeps); //convergence thresholds
    auto maxsweepsgroup = InputGroup(argv[1],"maxsweeps");
    auto maxsweepsvals = ReadIntVec(maxsweepsgroup,numeps); //max numbers of sweeps
    auto Noisegroup = InputGroup(argv[1],"Noise");
    auto Noisevals = ReadRealVec(Noisegroup,numeps); //noises (not necessary)
    auto minsweepsgroup = InputGroup(argv[1],"minsweeps");
    auto minsweepsvals = ReadIntVec(minsweepsgroup,numeps); //min numbers of sweeps
    auto MinDimgroup = InputGroup(argv[1],"MinDim");
    auto MinDimvals = ReadIntVec(MinDimgroup,numeps); //min bond dimensions
    
    // Final measurements
    auto storefinalmeas = parameters.getYesNo("savefinalmeas",true);
    // energies, bond dimensions, sweep numbers, computation times
    auto storewf = parameters.getYesNo("savewf",false); //MPS in ITensor format
    auto storewfMMA = parameters.getYesNo("savewfMMA",false); //Mathematica format
    auto storespcorr = parameters.getYesNo("savespcorr",true); //one-particle corr
    auto storennavg = parameters.getYesNo("savennavg",true); //<density-density>
    auto storeentropy = parameters.getYesNo("saveentropy",true); //entanglement
    auto storelocaldms = parameters.getYesNo("savelocaldms",true); //reduced density matrices (dms) for each segment
    
    // Intermediate measurements
    auto saveenergyhistory = parameters.getYesNo("saveenergyhistory",true); //energy and discontinuity after each DMRG step (local minimization)
    auto storelocalhistoryeverysweep = parameters.getYesNo("store_local_history",false); //dms after each sweep
    auto storelocalhistoryeveryeps = parameters.getYesNo("store_local_history_eps",false); //dms after each cycle
    
    // Product state in input file - string of integers (in which basis state)
    auto initialvec = InputGroup(argv[1],"initialvec");
    auto initialstate = ReadStringVec(initialvec,numsites);
    println("Got product state in input file");
    
    
//--------------------------------Create Site Set---------------------------------//
    
    
    auto store = SiteStore(numsites); //create sites storage
    
    auto numstates = parameters.getInt("numstates"); //number of basis states
    auto qngroup = InputGroup(argv[1],"qns"); //corresponding particle numbers
    qngroup.GotoGroup();
    std::vector<int> qns;
    for(int jj=1; jj<=numstates; ++jj)
    {
        int qn;
        qngroup.file() >> qn;
        qns.push_back(qn);
    }
    auto psiLrules = getRules(argv[1],"psiL"); //matrix elems of psi at left edge
    auto psiRrules = getRules(argv[1],"psiR"); //matrix elems of psi at right edge
    auto maxorder = parameters.getInt("maxorder"); //1 + max power of x in psi(x) or density(x)
    std::vector<int> nullqns(maxorder,0); //no qn associated with these powers
    auto psirules = getRulesVec(argv[1],"psi"); //array elems of psi(x) polynommial
    auto nrules = getRulesVec(argv[1],"n"); //array elems of density(x) polynommial
    
    for(int j=1; j<=numsites; ++j)
    {
        auto Hrules = getRules(argv[1],"H"+str(j)); //local Hamiltonians (H_c)
        auto siteobj = cSite(qns,psiLrules,psiRrules,Hrules,nullqns,psirules,nrules,{"SiteNumber",j}); //sites object
        store.set<itensor::cSite>(j,std::move(siteobj)); //move to site storage
    }
    // Create site set - which will be used to create MPS and MPO
    auto sites = BSiteSet(std::move(store));
    
    
//-------------------------------Make Hamiltonian---------------------------------//
    
    
     print("Making Hamiltonian...");
     auto ampo = AutoMPO(sites);
     for(int b=1; b<=numsites; ++b)
     {
         ampo += 1,"H",b; //"on-site" Hamiltonian
     }
     auto H = toMPO(ampo,{"Exact",true});
     println("done");
     
     
     print("Making penalty operator...");
     auto pampo = AutoMPO(sites);
     for(int b=1; b<=numsites; ++b)
     {
         pampo += 1,"psiRdagpsiR",b; //"on-site" terms
         pampo += 1,"psiLdagpsiL",b;
     }
     for(int b=1; b<numsites; ++b)
     {
         pampo += -1,"psiRdag",b,"psiL",b+1; //nearest-neighbor terms
         pampo += -1,"psiR",b,"psiLdag",b+1;
     }
     auto penaltyOp = toMPO(pampo,{"Exact",true}); //total discontinuity
     println("done");
    

//----------------------------Make/load initial state-----------------------------//
    
    
    MPS psi;
    auto cstate = InitState(sites); //constructor object for MPS
    
    if (loadfromfile) //load from a file, e.g., saved on a prior run
    {
        auto loadfilename = parameters.getString("loadfromfilename","default_psi");
        if (fileexists(loadfilename)) //check if file exists
        {
            print("Loading initial state from ",loadfilename);
            readFromFile(loadfilename,psi); //input MPS
            for (auto i=1; i<=numsites; ++i) //match tensor indices
            {
                auto snew = sites(i);
                auto sold = siteIndex(psi,i);
                auto newA = psi(i)*delta(snew,dag(sold));
                psi.set(i,newA);
            }
            goto initassigned; //state assigned
        }
    }
    
    print("Setting initial state from input file: "); //else, load product state
    for(int i=1; i<=numsites; ++i)
    {
        cstate.set(i,initialstate[i-1]);
        print(i,":",initialstate[i-1]);
        if (i<numsites) {print(",");}
    }
    psi = MPS(cstate); //state assigned
    
    initassigned: println("...done");
    psi.position(1); //right normalize
    
    
//--------------------------Initialize observer storage---------------------------//
    
    
    auto ehistname = rawname+"_ehist.m"; //name of energy history file
    auto localsweepname = rawname+"_localsweep.m"; //local density mats every sweep
    auto localepsname = rawname+"_localeps.m"; //local dms every eps cycle
    
    auto energy = inner(psi,H,psi); //initial energy
    auto discontinuity = inner(psi,penaltyOp,psi); //initial discontinuity
    Real energywithpenalty = NAN; //not a number - to be updated
    
    //
    // Write headers to storage files
    //
    if(saveenergyhistory) //store initial energy
    {
        print("Storing energy of initial state...");
        ofstream file;
        file.open(ehistname);
        file.precision(16); //16-digit precision
        file.setf(ios::fixed);
        file << "{{\"E\"->" << energy << ",\"disc\"->" << discontinuity;
        file << ",\"center\"->0,\"dir\"->0,"; //orthogonality center & sw direction
        file << "\"cputime\"->0,\"walltime\"->0,";
        file << "\"bonddim\"->0,\"eps\"->0}"; //local bond dimension & 1/penalty
        file.close();
        println("done");
    }
    
    if(storelocalhistoryeverysweep) //store dms for each segment after each sweep
    {
        print("Storing initial local dms for every sweep...");
        ofstream file;
        file.open(localsweepname);
        file.precision(16);
        file.setf(std::ios::fixed);
        file << "{";
        file.close();
        savelocal(localsweepname,psi,savethresh); //saveimg (add for imaginary part)
        println("done");
    }
    
    if (storelocalhistoryeveryeps) //store local dms after every eps cycle
    {
        ofstream file;
        print("Storing initial local dms for every eps...");
        file.open(localepsname);
        file.precision(16);
        file.setf(std::ios::fixed);
        file << "{";
        file.close();
        savelocal(localepsname,psi,savethresh); //saveimg //saves local dms
        println("done");
    }
    
    
//--------------------------Set initial sweep parameters--------------------------//
    
    
    print("Setting initial arguments...");
    auto epsargs = Args({
        "ShowEigs=",false, //details about singular value truncation
        "Truncate=",true, //truncates singular values using cutoff, mindim, maxdim
        "RespectDegenerate=",false, //degenerate subspaces all truncated or kept
        "PrintEigs=",true, //prints slowest-decaying eigenvalues after each sweep
        "Quiet=",false, //suppress most output except a short summary of each sweep
        "Silent=",false, //suppress all output and perform no measurements
        "StoreLocalHistory=",storelocalhistoryeverysweep, //dms after each sweep
        "LocalSweepName=",localsweepname, //file for storing dms after each sweep
        "SaveThresh=",savethresh, //don't store numbers smaller than this
        "saveenergyhistory=",saveenergyhistory, //energy after every step
        "ehistname=",ehistname, //file for storing energy history
        "DebugLevel",2, //prints extra info from eigensolver
        "UseSVD=",false //always use singular-value decomposition
    });//  Add others that are needed
    println("done");
    
    
//--------------------------------Find ground state-------------------------------//
    
    
    println("Starting eps sweeps"); //DMRG cycles with increasing penalties (1/eps)
    int sweep = 0;
    bool notconverged = true;
    std::vector<int> sweepnums; //number of sweeps for each cycle
    std::vector<double> finalens; //final energies after each cycle (E_c)
    std::vector<double> finaldiscs; //final discontinuities after each cycle
    std::vector<int> finalbonddims; //final bond dimensions after each cycle
    std::vector<double> walltimes; //wall times for each cycle
    std::vector<double> cputimes; //CPU times for each cycle
    
    std::ofstream lepsfile; //for storing local dms after each cycle
    if (storelocalhistoryeveryeps)
    {
        lepsfile.open(localepsname,std::ofstream::app); //open to append
        lepsfile.precision(16);
        lepsfile.setf(std::ios::fixed);
    }
    
    
    while (sweep<numeps and notconverged) //numeps = max number of sweeps
    {
        cpu_time eps_time; //measure duration of sweeps //alternatively, use clock()
        
        // Set sweep parameters
        epsargs.add("eps=",epsvals[sweep]); //value of eps (1/penalty)
        epsargs.add("MaxSweeps=",maxsweepsvals[sweep]); //max number of sweeps
        epsargs.add("MinSweeps=",minsweepsvals[sweep]); //min number of sweeps
        epsargs.add("SweepThresh=",threshvals[sweep]); //convergence threshold
        
        auto epssweeps = Sweeps(maxsweepsvals[sweep]); //sweeps object
        epssweeps.maxdim() = MaxDimvals[sweep]; //max bond dimension
        epssweeps.mindim() = MinDimvals[sweep]; //min bond dimension
        epssweeps.noise() = Noisevals[sweep]; //noise (usually set to 0)
        epssweeps.cutoff() = cutoffvals[sweep]; //singular-value cutoff
        epssweeps.niter() = MaxItervals[sweep]; //max eigensolver iterations
        
        // DMRG sweep for a given penalty
        energywithpenalty = cDMRGeps(psi,H,penaltyOp,epssweeps,epsargs); //sweeps
        discontinuity = epsargs.getReal("discontinuity"); //final discontinuity
        energy = epsargs.getReal("energy"); //final energy (E_c)
        
        int actualsw = epsargs.getInt("numsweeps"); //actual number of sweeps
        sweepnums.push_back(actualsw); //store as a list
        finalens.push_back(energy); //store energy
        finaldiscs.push_back(discontinuity); //store discontinuity
        finalbonddims.push_back(maxLinkDim(psi)); //store bond dimensions
        
        auto interval = eps_time.sincemark(); //cycle duration
        cputimes.push_back(interval.time); //store CPU time
        walltimes.push_back(interval.wall); //store wall time
        
        if (storelocalhistoryeveryeps) //store local dms after each eps
        {
            lepsfile << ",";
            savelocal(lepsfile,psi,savethresh); //saveimg
        }
        
        if (discontinuity<targetdiscontinuity) //stop if target is reached
        {
            notconverged=false;
            std::cout<<"Converged to discontinuity="<<discontinuity<<" with target="<<targetdiscontinuity<<"\n";
        }
        
        sweep+=1; //go to next eps (higher penalty)
        
    } //done with all eps sweeps
    
    
    if (storelocalhistoryeveryeps) //close local dms file for each cycle
    {
        lepsfile << "}";
        lepsfile.close();
    }
    
    if (notconverged) //print final discontinuity
    {
        std::cout<<"Did not converge: discontinuity="<<discontinuity<<" target="<<targetdiscontinuity<<"\n";
    }
    
    println("final energy = ",energy);
    println("bond dimension = ",maxLinkDim(psi),"\n");
    
    
//------------------------------Close observer storage----------------------------//
    
    
    //
    // Write tails to storage files
    //
    if(saveenergyhistory) //energy and discontinuity after each step
    {
        print("Closing energy observer file...");
        ofstream file;
        file.open(ehistname,std::ofstream::app);
        file << "}";
        file.close();
        println("done");
    }
    
    if(storelocalhistoryeverysweep) //local dms after each sweep
    {
        print("Closing localdm observer after each sweep...");
        ofstream file;
        file.open(localsweepname,std::ofstream::app);
        file << "}";
        file.close();
        println("done");
    }
    
    
//-----------------------------Save final measurements----------------------------//
    
    
    println("Final measurements");
    
    if (storespcorr) //save single-particle correlations (polynomial coefficients)
    {
        print("Saving single-particle correlation coefficients...");
        ofstream spcorrfile;
        spcorrfile.open(rawname+"_spcorrcoef.m");
        spcorrfile.precision(16);
        spcorrfile.setf(ios::fixed);
        savecorrvec(spcorrfile,psi,sites,"psidag","psi",savethresh); //saveimg
        spcorrfile.close();
        println("done");
    }
    
    if (storennavg) //save density-density correlations (polynomial coefficients)
    {
        print("Saving density-density average coefficients...");
        ofstream nnavgfile;
        nnavgfile.open(rawname+"_nnavgcoef.m");
        nnavgfile.precision(16);
        nnavgfile.setf(ios::fixed);
        savecorrvec(nnavgfile,psi,sites,"n","n",savethresh);
        nnavgfile.close();
        println("done");
    }
    
    if(storelocaldms) //save local density matrices
    {
        print("Saving local density matrices...");
        auto localdmfile=rawname+"_localdms.m";
        std::remove(localdmfile.c_str()); //erase file
        savelocal(localdmfile,psi,0.0); //saveimg //0.0 to compare small weights
        println("done");
    }
    
    if (storeentropy) //save von Neumann entanglement entropy for each bond
    {
        print("Saving von Neumann entropy...");
        ofstream entropyfile;
        entropyfile.open(rawname+"_entropy.m");
        entropyfile.precision(16);
        entropyfile.setf(ios::fixed);
        saveentropy(entropyfile,psi);
        entropyfile.close();
        println("done");
    }
    
    if (storewf) //save final state in ITensor format
    {
        print("Saving final MPS in ITensor format...");
        psi.position(1); //right normalize
        writeToFile(rawname+"_psi",psi); //Output MPS
        println("done");
    }
    
    if (storewfMMA) //save final MPS in Mathematica format (sparse array)
    {
        print("Saving final state for Mathematica...");
        ofstream psifile;
        psifile.open(rawname+"_psi_MMAformat.m");
        psifile.setf(ios::fixed);
        psifile.precision(16);
        psi.position(1); //right normalize
        savewfMMA(psifile,psi,savethresh); //saveimg
        psifile.close();
        println("done");
    }
    
    //
    // Computation time
    //
    auto overall_interval = overall_time.sincemark(); //total duration
    Real CPU_time = overall_interval.time; //total CPU time
    Real Wall_time = overall_interval.wall; //total wall time
    
    // Alternative way to measure total duration - used in benchmarking:
    // std::clock_t c_f = std::clock(); //final CPU time
    // auto t_f = std::chrono::high_resolution_clock::now(); //final wall time
    // Real CPU_time = 1.0 * (c_f-c_i) / CLOCKS_PER_SEC;
    // Real Wall_time = std::chrono::duration<double>(t_f-t_i).count();
    
    std::cout << "Total CPU time used: " << CPU_time << " s\n"
              << "Total wall time passed: " << Wall_time << " s\n";
    
    
    
    if (storefinalmeas) //save other measurements
    {
        print("Saving other measurements...");
        ofstream finalfile;
        finalfile.open(rawname+"_finalmeas.m");
        finalfile.precision(16);
        finalfile.setf(ios::fixed);
        finalfile << "{\"sweepnums\"->{"; //number of sweeps taken for each cycle
        for (int j=0; j<numeps-1; ++j) { finalfile << sweepnums[j] << ","; }
        finalfile << sweepnums[numeps-1] << "}";
        finalfile << ",\"energies\"->{"; //final energies after each cycle
        for (int j=0; j<numeps-1; ++j) { finalfile << finalens[j] << ","; }
        finalfile << finalens[numeps-1] << "}";
        finalfile << ",\"discontinuities\"->{"; //final discontinuities
        for (int j=0; j<numeps-1; ++j) { finalfile << finaldiscs[j] << ","; }
        finalfile << finaldiscs[numeps-1] << "}";
        finalfile << ",\"bonddims\"->{"; //final bond dimensions
        for (int j=0; j<numeps-1; ++j) { finalfile << finalbonddims[j] << ","; }
        finalfile << finalbonddims[numeps-1] << "}";
        finalfile << ",\"sweeptimesCPU\"->{"; //CPU times for each cycle
        for (int j=0; j<numeps-1; ++j) { finalfile << cputimes[j] << ","; }
        finalfile << cputimes[numeps-1] << "}";
        finalfile << ",\"sweeptimesWall\"->{"; //Wall times for each cycle
        for (int j=0; j<numeps-1; ++j) { finalfile << walltimes[j] << ","; }
        finalfile << walltimes[numeps-1] << "}";
        finalfile << ",\"Wall time\"->" << Wall_time; //total wall time
        finalfile << ",\"CPU time\"->" << CPU_time << "}"; //total CPU time
        finalfile.close();
        println("done");
    }
    
    
    return 0;
}
