//
//  cio.cpp
//  
//
//  Created by Erich Mueller on 6/13/19.
//  Edited by Shovan Dutta on 5/4/20.
//  Edited 7/28/21
//

#include "itensor/all.h"
#include "cio.h"


using namespace std;

namespace itensor {

    // Read sparse-array elements of the form <int, int, double>
    std::vector<inttriple> getRules(std::string filename,std::string groupname)
    {
        auto group = InputGroup(filename,groupname);
        group.GotoGroup();
        int numentries;  //number of nonzero elements
        group.file() >> numentries;
        std::vector<inttriple> rules;
        for(int jj=1; jj<=numentries; ++jj)
        {
            inttriple rule;
            group.file()>>get<0>(rule)>>get<1>(rule)>>get<2>(rule); //read elements
            rules.push_back(rule); //store elements
        }
        return rules;
    }
    
    // Read sparse-array elements of the form <int, int, int, double> [e.g. psi(x)]
    std::vector<quadruple> getRulesVec(std::string filename,std::string groupname)
    {
        auto group = InputGroup(filename,groupname);
        group.GotoGroup();
        int numentries;
        group.file() >> numentries;
        std::vector<quadruple> rules;
        for(int jj=1; jj<=numentries; ++jj)
        {
            quadruple rule;
            group.file()>>get<0>(rule)>>get<1>(rule)>>get<2>(rule)>>get<3>(rule);
            rules.push_back(rule);
        }
        return rules;
    }


    // Save two-site(segment) correlations <O^1_i O^2_j> for local ops O^1, O^2
    void savecorr(ofstream& file, MPS& psi, SiteSet& sites, string opname1, string opname2) //, bool saveimg //add this to save any imaginary parts
    {
        auto L = length(psi);
        file << "{";
        for(auto i=1; i<=L; ++i) //outer loop over segments
        {
            file << "{";
            auto leftop = op(sites,opname1,i);
            auto rightop = op(sites,opname2,i);
            psi.position(i); //set orthocenter at i
            auto psidag = prime(dag(psi),"Link");
            auto corr = prime(psidag(i),"Site")*leftop; // <O^1_i
            auto temp = prime(psi(i))*swapPrime(rightop,0,1); // O^2_i>
            auto val = elt(temp*corr); //eltC //change if complex
            file << val; //real(val)
            //if (saveimg) {file<<showpos<<imag(val)<<"I"<<noshowpos;}
            auto leftlink = leftLinkIndex(psi,i);
            corr *= prime(psi(i),leftlink); //contract with only O^1 at site i
            for(auto j=i+1; j<=L; ++j) //inner loop over segments
            {
                rightop = op(sites,opname2,j);
                corr *= psi(j);
                auto rightlink = commonIndex(corr,psi(j+1),"Link");
                auto corrval = prime(corr,rightlink)*rightop;
                corrval *= prime(psidag(j),"Site"); // <O^1_i O^2_j>
                val = elt(corrval); //eltC
                file<<","<<val; //real(val)
                //if (saveimg) {file<<showpos<<imag(val)<<"I"<<noshowpos;}
                if (j == L) break;
                corr *= psidag(j); //contract without O^2 at j
            }
            file << "}";
            if (i == L) break;
            file << ",";
        }
        file << "}";
    }

 
    // Save von Neumann entanglement entropy at all bonds (segment boundaries)
    void saveentropy(ofstream& file, MPS& psi)
    {
        auto L = length(psi);
        file << "{";
        for(auto b=1; b<L; ++b) //loop over bonds
        {
            psi.position(b); //ortho center at b
            auto l = leftLinkIndex(psi,b);
            auto s = siteIndex(psi,b);
            auto [U,S,V] = svd(psi(b),{l,s}); //singular-value decomposition
            auto u = commonIndex(U,S);
            Real SvN = 0.;
            for(auto n : range1(dim(u))) //number of singular values
            {
                auto Sn = elt(S,n,n);
                auto p = sqr(Sn); //square singular values
                if(p > 0.) SvN += -p*log(p); //add if nonzero
            }
            file << SvN;
            if (b == L-1) break;
            file << ",";
        }
        file << "}";
    }

    
    // Save correlations <O^1_i O_2_j> for operators with an additional index representing polynomial coefficients - e.g., psi(x) and density(x)
    void savecorrvec(ofstream& file, MPS& psi, SiteSet& sites, string opname1, string opname2, Real savethresh) //, bool saveimg //add if complex
    {
        auto L = length(psi);
        file << "{";
        for(auto i=1; i<=L; ++i) //outer loop over segments
        {
            file << "{";
            auto leftop = op(sites,opname1,i);
            auto rightop = swapPrime(op(sites,opname2,i),0,1);
            psi.position(i); //set ortho center at i
            auto psidag = prime(dag(psi),"Link");
            auto corr = prime(psidag(i),"Site")*leftop; // <O^1_i
            auto temp = prime(psi(i))*rightop; // O^2_i>
            auto array = temp*corr; // <O^1_i O^2_i> - matrix of coefficients
            auto leftind = findIndex(leftop,"Vector");
            auto rightind = findIndex(rightop,"Vector");
            savearray(file,permute(array,{leftind,rightind}),savethresh); //saveimg //save coefficient array to file
            auto leftlink = leftLinkIndex(psi,i);
            corr *= prime(psi(i),leftlink); //contract with only O^1 at site i
            for(auto j=i+1; j<=L; ++j) //inner loop over segments
            {
                rightop = prime(op(sites,opname2,j),"Vector");
                corr *= psi(j);
                auto rightlink = commonIndex(corr,psi(j+1),"Link");
                auto corrval = prime(corr,rightlink)*rightop;
                corrval *= prime(psidag(j),"Site"); // <O^1_i O^2_j> - matrix
                rightind = findIndex(rightop,"Vector");
                file<<",";
                savearray(file,permute(corrval,{leftind,rightind}),savethresh); //saveimg //save to file
                if (j == L) break;
                corr *= psidag(j); //contract without O^2 at j
            }
            file << "}";
            if (i == L) break;
            file << ",";
        }
        file << "}";
    }


    // Save MPS as a sparse array readable by Mathematica
    void savewfMMA(ofstream& psifile, MPS& psi, Real savethresh) //, bool saveimg //add if complex
    {
        auto L = length(psi);
        psifile << "{";
        for(auto j=1; j<=L; ++j) //loop over segments
        {
            IndexSet is; //ordered indices
            if (j==1)
            {
                is = {siteIndex(psi,1),rightLinkIndex(psi,1)}; //first site
            }
            else if (j==L)
            {
                is = {siteIndex(psi,L),leftLinkIndex(psi,L)}; //last site
            }
            else
            {
                is = {siteIndex(psi,j),leftLinkIndex(psi,j),rightLinkIndex(psi,j)}; //other sites
            }
            
            savearray(psifile,permute(psi(j),is),savethresh); //saveimg //save local array to file
            
            if(j == L) break;
            psifile<<",";
        }
        psifile<<"}";
    }


    // Save a (sparse) array in a Mathematica format
    void savearray(ofstream& file, ITensor dm, Real savethresh) //, bool saveimg
    {
        auto allinds = dm.inds();
        auto numinds = order(allinds);
        auto integerindices = std::vector<int>(numinds); //indices we print out
        auto indlens = std::vector<int>(numinds); //number of values each index takes on
        for (int i=0; i<numinds; ++i)
        {
            integerindices[i] = 1;  //Set default
            indlens[i] = allinds[i].size(); //Get length of each index value
        }
        
        bool happy = true;
        bool empty = true;
        
        file<<"SparseArray[{";
        
        while(happy)
        {
            auto val=dm.elt(integerindices); //eltC //change if complex
            
            if (abs(val)>savethresh) //only store if greater than a threshold
            {
                if(not(empty)) // not first element
                {
                    file<<","; //separator
                }
                file<<"{";
                for (int j=0;j<numinds;j++)
                {
                    file<<integerindices[j];
                    if (j<numinds-1)
                    {
                        file<<",";
                    }
                }
                file<<"}->"<<val; //real(val) //change if complex
                //if (saveimg) {file<<showpos<<imag(val)<<"I"<<noshowpos;}
                empty=false;
            }
            
            
            // increment indices
            int level=0;
            while(integerindices[level]==indlens[level]) // enter loop if counter is at max
            {
                level+=1;
                if (level>= numinds)
                {
                    break;
                }
            }
            if (level>=numinds) // gone through all of the indices
            {
                break;
            }
            integerindices[level]+=1;
            for(int j=0;j<level;j++)
            {
                integerindices[j]=1;
            }
        } // done adding elements
        
        file << "},{";
        for (int i=0;i<numinds;i++)
        {
            file<<indlens[i];
            if (i<numinds-1)
            {
                file<<",";
            }
        }
        file <<"}]";
    }

    
    // Save reduced density matrices for each segment
    void savelocal(std::string const& filename, MPS& psi, Real savethresh) //, bool saveimg //add if complex
    {
        ofstream file;
        file.open(filename,std::ofstream::app);  //open file for appending
        file.precision(16); //16-digit precision
        file.setf(ios::fixed);
        savelocal(file,psi,savethresh); //saveimg
        file.close();
    }
    
    
    void savelocal(std::ofstream& localfile, MPS& psi, Real savethresh) //, bool saveimg
    {
        auto numsites=length(psi);
        localfile << "{";
        for(auto j=1;j<=numsites;j++) //loop over segments
        {
            psi.position(j); //ortho center at j
            auto A=psi.A(j); //local tensor
            auto dA=dag(prime(psi(j),"Site"));
            auto dm=dA*A; //density matrix
            savearray(localfile,dm,savethresh); //saveimg //save to file
            if(j<numsites)
            {
                localfile<<",";
            }
        }
        localfile <<"}";
    }


} //end namespace
