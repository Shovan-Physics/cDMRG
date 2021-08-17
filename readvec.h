//
//  readvec.h
//
//
//  Created by Erich Mueller on 3/6/20.
//

#ifndef readvec_h
#define readvec_h

#include "itensor/util/input.h"

namespace itensor
{

    std::vector<Real> ReadRealVec(InputGroup &group, int len);
    std::vector<std::string> ReadStringVec(InputGroup &group, int len);
    std::vector<int> ReadIntVec(InputGroup &group, int len);

}

#endif /* readvec_h */
