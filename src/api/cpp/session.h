//
// Created by Eric Irrgang on 11/29/17.
//

#ifndef GROMACS_SESSION_H
#define GROMACS_SESSION_H

namespace gmxapi
{

// forward declaration
class Status;

class Session
{
    public:
        Status run();
};

}

#endif //GROMACS_SESSION_H
