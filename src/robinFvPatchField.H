/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#ifndef robinFvPatchField_H
#define robinFvPatchField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class robinFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template <class Type>
class robinFvPatchField
    : public mixedFvPatchField<Type>
{
protected:
    // Private data

    // exchange coeff (dT/dt = beta_(T_inf - T))
    scalar beta_;

public:
    //- Runtime type information
    TypeName("robin");

    // Constructors

    //- Construct from patch and internal field
    robinFvPatchField(
        const fvPatch &,
        const DimensionedField<Type, volMesh> &);

    //- Construct from patch, internal field and dictionary
    robinFvPatchField(
        const fvPatch &,
        const DimensionedField<Type, volMesh> &,
        const dictionary &);

    //- Construct by mapping given robinFvPatchField
    //  onto a new patch
    robinFvPatchField(
        const robinFvPatchField<Type> &,
        const fvPatch &,
        const DimensionedField<Type, volMesh> &,
        const fvPatchFieldMapper &);

    //- Construct as copy
    robinFvPatchField(
        const robinFvPatchField &);

    //- Construct and return a clone
    virtual tmp<fvPatchField<Type>> clone() const
    {
        return tmp<fvPatchField<Type>>(
            new robinFvPatchField<Type>(*this));
    }

    //- Construct as copy setting internal field reference
    robinFvPatchField(
        const robinFvPatchField &,
        const DimensionedField<Type, volMesh> &);

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchField<Type>> clone(
        const DimensionedField<Type, volMesh> &iF) const
    {
        return tmp<fvPatchField<Type>>(
            new robinFvPatchField<Type>(*this, iF));
    }

    // Member functions

    // Access

    // Evaluation functions

    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();

    //- Write
    virtual void write(Ostream &) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "robinFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
