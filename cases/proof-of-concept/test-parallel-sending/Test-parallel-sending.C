/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    Test-parallel-sending

Description
    Test the procedure of every processor sending to everyone else.
    This will be needed when OpenFOAM and Lammps partitions are different.

\*---------------------------------------------------------------------------*/

#include "List.H"
#include "mapDistribute.H"
#include "argList.H"
#include "Time.H"
#include "IPstream.H"
#include "OPstream.H"
#include "vector.H"
#include "IOstreams.H"
#include "Random.H"
#include "Tuple2.H"
#include "PstreamBuffers.H"
#include "vectorList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"

    // get the number of processors
    label nproc=Pstream::nProcs();

    vector localV = vector::zero;

    label myrank = Pstream::myProcNo();

    // initialize lists for sending and receiving with mask initial values
    vectorList myNos(nproc, vector(555, 555, 555));
    vectorList receivedNos(nproc, vector(666, 666, 666));

    Info << "Total # of procs: " << nproc << endl;
    

    for (int i=0; i<nproc; i++)
    {
        scalar s = i + myrank*10;
        myNos[i] = vector(s, 100+s, 200+s);
    }

    Pout << "My initial list: " << myNos << endl;


    PstreamBuffers pBufs(Pstream::nonBlocking);

    // send to everyone else with UOPstream (output)
    for(int i=0; i<nproc; i++)
    {
        if (i != myrank)
        {
            UOPstream toEveryoneElse(i, pBufs);
            toEveryoneElse << myNos[i];
        }
    }

    // finish sending and block threads here
    pBufs.finishedSends();

    // receive from everyone else; copy local data residing on my proc
    for(int i=0; i<nproc; i++)
    {
        if (i != myrank)
        {
            UIPstream fromEveryoneElse(i, pBufs);
            receivedNos[i] = vector(fromEveryoneElse);
        }
        else
        {
            receivedNos[i] = myNos[i];
        }
    }


    /* Another way of achieving the same goal */ 
    
    //    Pstream::resetRequests(0);
    // // Send the i^th number to the i^th processor
    // // omit if sending to myself, i.e., i = myrank
    // for (int i=0; i<nproc; i++)
    // {
    //     if (i != myrank)
    //     {
    //         OPstream vectorStream(Pstream::blocking, i);
    //         vectorStream << myNos[i];
    //     }
    // }

    // vectorStream.finishedSends();

    // // Create buffers for recieiving messages
    // for (int i=0; i<nproc; i++)
    // {
    //     if (i != myrank)
    //     {
    //         IPstream vStream(Pstream::blocking, i);
    //         vStream >> receivedNos[i];
    //     }
    //     //        else
    //     //        {
    //         //receivedNos[i] = myNos[i];
    //     //        }
    // }

    // UPstream::waitRequests();
    
    Pout << "My received list: " << receivedNos << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
