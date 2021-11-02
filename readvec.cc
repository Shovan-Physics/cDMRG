//
//  readvec.cpp
//
//
//  Created by Erich Mueller on 3/6/20.
//

#include <readvec.h>
#include "itensor/util/input.h"

namespace itensor
{

    std::vector<Real> ReadRealVec(InputGroup &group, int len) //vector of real numbers
    {
        std::vector<Real> result;
        group.GotoGroup();
        for (int j = 1; j <= len; j++)
        {
            Real val;
            group.file() >> val;
            result.push_back(val);
        }
        return result;
    }

    std::vector<std::string> ReadStringVec(InputGroup &group, int len) //strings
    {
        std::vector<std::string> result;
        group.GotoGroup();
        for (int j = 1; j <= len; j++)
        {
            std::string val;
            group.file() >> val;
            result.push_back(val);
        }
        return result;
    }

    std::vector<int> ReadIntVec(InputGroup &group, int len) //integers
    {
        std::vector<int> result;
        group.GotoGroup();
        for (int j = 1; j <= len; j++)
        {
            int val;
            group.file() >> val;
            result.push_back(val);
        }
        return result;
    }

} // End Namespace
