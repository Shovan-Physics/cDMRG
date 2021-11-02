//
//  cio.hpp
//
//
//  Created by Erich Mueller on 6/13/19.
//  Edited by Shovan Dutta on 5/4/20.
//  Edited 7/28/21
//

#ifndef cio_hpp
#define cio_hpp

#include <stdio.h>
#include "itensor/all.h"

namespace itensor
{

    // New structure for array elements - same as in site set (cSite.h)
    using inttriple = std::tuple<int, int, double>;
    using quadruple = std::tuple<int, int, int, double>;

    std::vector<inttriple> getRules(std::string filename, std::string groupname);
    std::vector<quadruple> getRulesVec(std::string filename, std::string groupname);

    void savecorr(std::ofstream &file, MPS &psi, SiteSet &sites, string opname1, string opname2); //, bool saveimg=true); //add if complex
    void saveentropy(std::ofstream &file, MPS &psi);
    void savecorrvec(std::ofstream &file, MPS &psi, SiteSet &sites, string opname1, string opname2, Real savethresh); //, bool saveimg=true);
    void savewfMMA(std::ofstream &psifile, MPS &psi, Real savethresh);                                                //, bool saveimg=true);
    void savearray(std::ofstream &file, ITensor dm, Real savethresh);                                                 //, bool saveimg=true);
    void savelocal(std::string const &filename, MPS &psi, Real savethresh);                                           //, bool saveimg=true);
    void savelocal(std::ofstream &file, MPS &psi, Real savethresh);                                                   //, bool saveimg=true);

}
#endif /* cio_hpp */
