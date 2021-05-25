/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    dispersionEvaluationFoam

Description
  Solver for the assessment of dispersison tensor of a porous structure
  Solves the closure problem by Carbonell and Whitaker 1983 to evaluate
  dispersion tensor from pore-scale simulations on a cyclic REV given a
  local velocity profile U.

Author
  05-01-2021 - CS - upgrade to OpenFOAM 7
	31-12-2013 - Cyprien Soulaine (CS) - cyprien.soulaine@gmail.com

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"

//#include "BboundaryCondition/BboundaryFvPatchVectorField.H"
#include "fBoundaryCondition/fBoundaryConditionFvPatchVectorField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "CourantNo.H"


    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            fvVectorMatrix BEqn
            (
                (U-Umoy)
        //      +  fvm::ddt(B)
              + fvm::div(phi, B)
              - fvm::laplacian(Dmol, B)
            );

            BEqn.setReference(BRefCell, BRefValue);
            BEqn.relax();
            BEqn.solve();

        }
        //b = B-B.weightedAverage(mesh.V());
        b = B -fvc::domainIntegrate(B)/sum(mesh.V());
//        volVectorField b ("b", B-B.weightedAverage(mesh.V()));

        runTime.write();
        if(runTime.outputTime())
        {
            b.write();
        }
    }


    dimensionedTensor D ("D", dimensionSet(0,0,0,0,0,0,0), tensor::zero);

//    volVectorField b ("b", B-B.weightedAverage(mesh.V()));
//    b.write();

    volTensorField Utildab ("Utildab",((U-Umoy)*b));

    volTensorField gradb ("gradb", fvc::grad(b));

    D = tensor(1,0,0,0,1,0,0,0,1)*Dmol/Dmol - Utildab.weightedAverage(mesh.V())/Dmol + gradb.weightedAverage(mesh.V());


    Info<< "D / Dmol = " << D.value() << nl << endl;


    Info<< "Dzz / Dmol = " << D.component(tensor::ZZ).value() << nl << endl;


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
