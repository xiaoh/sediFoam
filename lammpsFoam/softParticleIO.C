/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is based on OpenFOAM.

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

\*---------------------------------------------------------------------------*/

#include "softParticleCloud.H"
#include "IOstreams.H"
#include "softParticle.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::softParticle::softParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    particle(mesh, is, readFields),
    d_(0.0),
    mass_(0.0),
    U_(vector::zero),
    moveU_(vector::zero),
    ensembleU_(vector::zero),
    positionOld_(vector::zero),
    UOld_(vector::zero),
    tag_(0),
    lmpCpuId_(0),
    type_(0),
    density_(0.0)
{
    // Pout<< "Creating a particle from Istream.." << endl;

    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            d_ = readScalar(is);
            is >> U_;
            is >> moveU_;
            is >> ensembleU_;
            is >> tag_;
            is >> lmpCpuId_;
            is >> type_;
            density_ = readScalar(is);
        }
        else
        {
            is.read
            (
               reinterpret_cast<char*>(&d_),
               sizeof(mass_) + sizeof(d_)  + sizeof(positionOld_) + sizeof(UOld_)
             + sizeof(U_) + sizeof(moveU_) + sizeof(ensembleU_) + sizeof(type_)
             + sizeof(density_) + sizeof(tag_)  + sizeof(lmpCpuId_)
            );
          }
    }

    // calculateDerived();
    if (debug)
    {
        Pout<< "creating a softParticle." << endl;
        Pout<< "tag is: " << tag_ << endl;
        Pout<< "lmpCpuId is: " << lmpCpuId_ << endl;
        Pout<< "type is: " << type_ << endl;

        Pout<< "d is: " << d_ << endl;
        Pout<< "U is: " << U_ << endl;
        Pout<< "moveU is: " << moveU_ << endl;
        Pout<< "ensembleU is: " << ensembleU_ << endl;
        Pout<< "density is: " << density_ << endl;
        Pout<< "mass is: " << mass_ << endl;
        Pout<< "positionOld is: " << positionOld_ << endl;
        Pout<< "UOld is: " << UOld_ << endl;
        Pout<< "position is: " << position() << endl;
    }

    // Check state of Istream
    is.check("softParticle::softParticle(Istream&)");
}


void Foam::softParticle::readFields(Cloud<softParticle>& c)
{
    if (!c.size())
    {
        return;
    }

    particle::readFields(c);

    IOField<scalar> d(c.fieldIOobject("d", IOobject::MUST_READ));
    c.checkFieldIOobject(c, d);

    IOField<scalar> density(c.fieldIOobject("density", IOobject::MUST_READ));
    c.checkFieldIOobject(c, density);

    IOField<scalar> tag(c.fieldIOobject("tag", IOobject::MUST_READ));
    c.checkFieldIOobject(c, tag);

    IOField<scalar> lmpCpuId(c.fieldIOobject("lmpCpuId", IOobject::MUST_READ));
    c.checkFieldIOobject(c, lmpCpuId);

    IOField<scalar> type(c.fieldIOobject("type", IOobject::MUST_READ));
    c.checkFieldIOobject(c, type);

    IOField<vector> U(c.fieldIOobject("U", IOobject::MUST_READ));
    c.checkFieldIOobject(c, U);

    label i = 0;
    forAllIter(Cloud<softParticle>, c, iter)
    {
        softParticle& p = iter();

        p.d_ = d[i];
        p.density_ = density[i];
        p.tag_ = tag[i];
        p.lmpCpuId_ = lmpCpuId[i];
        p.type_ = type[i];
        p.U_ = U[i];
        i++;
    }
}


void softParticle::writeFields(const Cloud<softParticle>& c)
{
    particle::writeFields(c);

    label np = c.size();

    IOField<scalar> d(c.fieldIOobject("d", IOobject::NO_READ), np);
    IOField<scalar> density(c.fieldIOobject("density", IOobject::NO_READ), np);
    IOField<label> tag(c.fieldIOobject("tag", IOobject::NO_READ), np);
    IOField<label> lmpCpuId(c.fieldIOobject("lmpCpuId", IOobject::NO_READ), np);
    IOField<label> type(c.fieldIOobject("type", IOobject::NO_READ), np);
    IOField<vector> U(c.fieldIOobject("U", IOobject::NO_READ), np);
    IOField<vector> ensembleU
    (
        c.fieldIOobject("ensembleU", IOobject::NO_READ),
        np
    );

    label i = 0;
    forAllConstIter(Cloud<softParticle>, c, iter)
    {
        const softParticle& p = iter();

        d[i] = p.d_;
        density[i] = p.density_;
        tag[i] = p.tag_;
        lmpCpuId[i] = p.lmpCpuId_;
        type[i] = p.type_;
        U[i] = p.U_;
        ensembleU[i] = p.ensembleU_;
        i++;
    }

    d.write();
    tag.write();
    lmpCpuId.write();
    type.write();
    U.write();
    ensembleU.write();
}


void softParticle::writeFields(const Cloud<softParticle>& c, const label np)
{
    IOField<vector> positions
    (
        c.fieldIOobject("positions", IOobject::NO_READ),
        np
    );
    IOField<scalar> d(c.fieldIOobject("d", IOobject::NO_READ), np);
    IOField<label> tag(c.fieldIOobject("tag", IOobject::NO_READ), np);
    IOField<label> lmpCpuId(c.fieldIOobject("lmpCpuId", IOobject::NO_READ), np);
    IOField<label> type(c.fieldIOobject("type", IOobject::NO_READ), np);
    IOField<vector> U(c.fieldIOobject("U", IOobject::NO_READ), np);

    IOField<vector> ensembleU
    (
        c.fieldIOobject("ensembleU", IOobject::NO_READ),
        np
    );

    IOField<vector> positionsOld
    (
        c.fieldIOobject("positionsOld", IOobject::NO_READ),
        np
    );

    IOField<vector> UOld
    (
        c.fieldIOobject("UOld", IOobject::NO_READ),
        np
    );

    IOField<label> origProc
    (
        c.fieldIOobject("origProcId", IOobject::NO_READ),
        np
    );

    IOField<label> origId(c.fieldIOobject("origId", IOobject::NO_READ), np);

    label i = 0;
    forAllConstIter(Cloud<softParticle>, c, iter)
    {
        const softParticle& p = iter();

        if (p.cell() < 0) continue;

        positions[i] = p.position();

        d[i] = p.d_;
        tag[i] = p.tag_;
        lmpCpuId[i] = p.lmpCpuId_;
        type[i] = p.type_;
        U[i] = p.U_;
        ensembleU[i] = p.ensembleU_;
        positionsOld[i] = p.positionOld_;
        UOld[i] = p.UOld_;

        origProc[i] = p.origProc();
        origId[i] = p.origId();
        i++;

}

    positions.write();
    positionsOld.write();
    d.write();
    tag.write();
    lmpCpuId.write();
    type.write();
    ensembleU.write();
    origProc.write();
    origId.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const softParticle& p)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.d_
            << token::SPACE << p.tag_
            << token::SPACE << p.lmpCpuId_
            << token::SPACE << p.type_
            << token::SPACE << p.U_
            << token::SPACE << p.moveU_
            << token::SPACE << p.ensembleU_
            << token::SPACE << p.mass_
            << token::SPACE << p.density_
            << token::SPACE << p.positionOld_
            << token::SPACE << p.UOld_;
    }
    else
    {
        os  << static_cast<const particle&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.d_),
            sizeof(p.d_) + sizeof(p.tag_) + sizeof(p.lmpCpuId_) + sizeof(p.type_) + sizeof(p.U_)
          + sizeof(p.moveU_) + sizeof(p.ensembleU_) + sizeof(p.mass_)
          + sizeof(p.density_) + sizeof(p.positionOld_) + sizeof(p.UOld_)
        );
    }

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const softParticle&)");

    return os;
}


// ************************************************************************* //

