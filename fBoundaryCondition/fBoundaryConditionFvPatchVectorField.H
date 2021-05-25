/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Class
    Foam::fBoundaryConditionFvPatchScalarField

Group
    grpGenericBoundaryConditions grpWallBoundaryConditions

Description
    This boundary condition applies at the fluid-solid boundary for solving the
    closure problems of Carbonell and Whitaker 1983.

    It solves:
        n.grad(B)=-n

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            closureBoundaryConditionForB;
        value		uniform (0 0 0);
    }
    \endverbatim

SourceFiles
    fBoundaryConditionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fBoundaryConditionFvPatchVectorFields_H
#define fBoundaryConditionFvPatchVectorFields_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class fBoundaryConditionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fBoundaryConditionFvPatchVectorField
:
    public fixedGradientFvPatchVectorField
{
    // Private data


public:

    //- Runtime type information
    TypeName("closureBoundaryConditionForB");


    // Constructors

        //- Construct from patch and internal field
        fBoundaryConditionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fBoundaryConditionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  fBoundaryConditionFvPatchScalarField onto a new patch
        fBoundaryConditionFvPatchVectorField
        (
            const fBoundaryConditionFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fBoundaryConditionFvPatchVectorField
        (
            const fBoundaryConditionFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new fBoundaryConditionFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fBoundaryConditionFvPatchVectorField
        (
            const fBoundaryConditionFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new fBoundaryConditionFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
