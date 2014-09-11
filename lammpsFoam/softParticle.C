/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*----------------------------------------------------------------------------*/

#include "softParticle.H"
#include "IOstreams.H"
#include "softParticleCloud.H"
#include "wallPolyPatch.H"
#include "processorPolyPatch.H"

namespace Foam {

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
softParticle::softParticle
(
    const polyMesh& mesh,
    const vector& position,
    const label celli,
    const scalar& d,
    const vector& U,
    const scalar& rhos,
    const label& tag,
    const label& lmpCpuId,
    const label& type
)
:
    particle(mesh,position,celli),
    d_(d),
    mass_(0.),
    U_(U),
    moveU_(vector::zero),
    ensembleU_(vector::zero),
    positionOld_(vector::zero),
    tag_(tag),
    lmpCpuId_(lmpCpuId),
    type_(type),
    density_(rhos)
{
    calculateDerived();
}

// Calculating some properties of the created particles.
void softParticle::calculateDerived()
{
    mass_=  density_*4./3.*constant::mathematical::pi*d_*d_*d_/8.;
    positionOld_ = position_;
    // Pout<< "creating a softParticle." << endl;
    if (debug)
    {
        Pout<< "tag is: " << tag_ << endl;
        Pout<< "type is: " << type_ << endl;
        Pout<< "d is: " << d_ << endl;
        Pout<< "U is: " << U_ << endl;
        Pout<< "density is: " << density_ << endl;
        Pout<< "mass is: " << mass_ << endl;
        Pout<< "positionOld is: " << positionOld_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

softParticle::~softParticle()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Move a softParticle from one position to another.
bool Foam::softParticle::move
(
    trackingData& td,
    const scalar trackTime
)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    const polyBoundaryMesh& pbMesh = mesh_.boundaryMesh();

    scalar tEnd = (1.0 - stepFraction())*trackTime;
    scalar dtMax = tEnd;

    while (td.keepParticle && !td.switchProcessor && tEnd > ROOTVSMALL)
    {
        if (debug)
        {
            Pout<< "Time = " << mesh_.time().timeName()
                << " trackTime = " << trackTime
                << " steptFraction() = " << stepFraction() << endl;
            Pout<< "Before trackToFace, particle tag is: " << ptag()
                << ". Particle position is: " << position() << endl;
        }

        scalar dt = min(dtMax, tEnd);

        dt *= trackToFace(position() + dt*moveU_, td);

        tEnd -= dt;
        stepFraction() = 1.0 - tEnd/trackTime;

        if (debug)
        {
            Pout<< "After trackToFace, particle tag is: " << ptag()
                << ". Particle position is: " << position() << endl;
        }

        if (onBoundary() && td.keepParticle)
        {
            if (isA<processorPolyPatch>(pbMesh[patch(face())]))
            {
                td.switchProcessor = true;
                // Pout<< "Cross the processor boundary.." << endl;
            }
        }
    }

    return td.keepParticle;
}


bool Foam::softParticle::hitPatch
(
    const polyPatch&,
    trackingData& td,
    const label,
    const scalar,
    const tetIndices&
)
{
    return false;
}


void Foam::softParticle::hitProcessorPatch
(
    const processorPolyPatch&,
    trackingData& td
)
{
    td.switchProcessor = true;
}


void Foam::softParticle::hitPatch
(
    const polyPatch&,
    trackingData& td
)
{
    td.keepParticle = false;
}

void Foam::softParticle::transformProperties (const tensor& T)
{
    particle::transformProperties(T);
    moveU_ = transform(T, moveU_);
    // Pout<< "hitting cyclic patch T" << endl;
}


void Foam::softParticle::transformProperties(const vector& separation)
{
    particle::transformProperties(separation);
    // Pout<< "hitting cyclic patch separation" << endl;
}



} // namespace Foam


// ************************************************************************* //
