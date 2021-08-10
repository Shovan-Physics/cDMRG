//
//  localmpoproj.hpp
//  
//
//  Created by Erich Mueller on 6/25/19.
//  Edited by Shovan Dutta on 7/28/21
//

#ifndef localmpoproj_hpp
#define localmpoproj_hpp

#include <stdio.h>
#include "itensor/mps/localmpo.h"

namespace itensor {
    
    class LocalMPOproj
    {
    private:
        MPO const* Op_ = nullptr; // Hamiltonian (H_c)
        MPO const* Proj_ = nullptr; // Discontinuity op - positive semidefinite
        LocalMPO lmpoH_;
        LocalMPO lmpoProj_;
        
    public:
        Real weight_ = 1; // How much to weight penalty
        
    LocalMPOproj() { }
        
    LocalMPOproj(MPO const& Op,
                 MPO const& Proj,
                 Args const& args = Args::global());
        
        //
        // Sparse matrix methods
        //
        
        void
        product(ITensor const& phi,
                ITensor& phip) const;
        
        Real
        expect(ITensor const& phi) const // Energy with penalty
        {
            return lmpoH_.expect(phi)+weight_*lmpoProj_.expect(phi);
        }
        
        Real
        Hexpect(ITensor const& phi) const // Energy (E_c)
        {
            return lmpoH_.expect(phi);
        }
        
        Real
        Projexpect(ITensor const& phi) const // Discontinuity measure
        {
            return lmpoProj_.expect(phi);
        }
        
        ITensor
        deltaRho(ITensor const& AA,
                 ITensor const& comb,
                 Direction dir) const
        {
            return lmpoH_.deltaRho(AA,comb,dir) + lmpoProj_.deltaRho(AA,comb,dir);
        }
        // note no weight factor!!
        // That was a design decision
        // can re-evaluate
        
        void
        position(int b, MPS const& psi);
        
        size_t
        size() const { return lmpoH_.size(); }
        
        explicit
        operator bool() const { return (bool(Op_) and bool(Proj_)); }
        
        Real
        weight() const { return weight_; }
        void
        weight(Real val) { weight_ = val; } // Set weight
        
    }; // end templates

    
    // Constructors
    
    inline LocalMPOproj::
    LocalMPOproj(MPO const& Op,
                 MPO const& Proj,
                 Args const& args)
    : Op_(&Op),Proj_(&Proj),weight_(args.getReal("Weight",1))
    {
        lmpoH_=LocalMPO(Op);
        lmpoProj_=LocalMPO(Proj);
    }

    void inline LocalMPOproj::
    product(ITensor const& phi,
            ITensor & phip) const
    {
        lmpoH_.product(phi,phip);
        ITensor other;
        lmpoProj_.product(phi,other);
        phip+=weight_*other;
    }
    
    void inline LocalMPOproj::
    position(int b, const MPS& psi)
    {
        lmpoH_.position(b,psi);
        lmpoProj_.position(b,psi);
    }


}
#endif /* localmpoproj_hpp */
