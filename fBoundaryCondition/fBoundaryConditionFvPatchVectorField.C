/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "fBoundaryConditionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //



Foam::fBoundaryConditionFvPatchVectorField::
fBoundaryConditionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF)
{}


Foam::fBoundaryConditionFvPatchVectorField::
fBoundaryConditionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF)
{
        fvPatchField<vector>::operator=(patchInternalField());
        gradient() = vector::zero;
}



Foam::fBoundaryConditionFvPatchVectorField::
fBoundaryConditionFvPatchVectorField
(
    const fBoundaryConditionFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::fBoundaryConditionFvPatchVectorField::
fBoundaryConditionFvPatchVectorField
(
    const fBoundaryConditionFvPatchVectorField& ptf
)
:
    fixedGradientFvPatchVectorField(ptf)
{}


Foam::fBoundaryConditionFvPatchVectorField::
fBoundaryConditionFvPatchVectorField
(
    const fBoundaryConditionFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fBoundaryConditionFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }


	gradient() = -(patch().nf());


    fixedGradientFvPatchVectorField::updateCoeffs();
}


void Foam::fBoundaryConditionFvPatchVectorField::write(Ostream& os) const
{
    fixedGradientFvPatchVectorField::write(os);
    writeEntry(os,"value",*this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        fBoundaryConditionFvPatchVectorField
    );
}

// ************************************************************************* //
