/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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

\*---------------------------------------------------------------------------*/

#include "linGeomVelocitySolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(linGeomVelocitySolid, 0);
addToRunTimeSelectionTable(physicsModel, linGeomVelocitySolid, solid);
addToRunTimeSelectionTable(solidModel, linGeomVelocitySolid, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linGeomVelocitySolid::linGeomVelocitySolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_),
    gradU_(fvc::grad(U()))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool linGeomVelocitySolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    // Mesh update loop
    do
    {
        int iCorr = 0;
        lduSolverPerformance solverPerfU;
        blockLduMatrix::debug = 0;

        Info<< "Solving the momentum equation for U" << endl;


        // Stablisation viscosity
        // const dimensionedScalar eta_ =
        //    solidProperties().lookupOrDefault<dimensionedScalar>
        //    (
        //        "numericalViscosity",
        //        dimensionedScalar("eta", dimless/dimTime, 0.0)
        //    );

        // Momentum equation loop
        do
        {
            // Store fields for under-relaxation and residual calculation
            //D().storePrevIter();


            U().storePrevIter();


            // Linear momentum equation total displacement form
//            fvVectorMatrix DEqn
//            (
//                rho()*fvm::d2dt2(D())
//                //+ eta_*rho_*fvm::ddt(D())
//             ==
//            fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
//              - fvc::laplacian(impKf_, D(), "laplacian(DD,D)")
//            + fvc::div(sigma(), "div(sigma)")
//              + rho()*g()
////              + mechanical().RhieChowCorrection(D(), gradD())
//            );

            fvVectorMatrix UEqn
            (
                rho()*fvm::ddt(U())
             == fvm::laplacian(runTime().deltaT()*impKf_, U())
              - fvc::laplacian(runTime().deltaT()*impKf_, U())
              + fvc::div(sigma())
              + rho()*g()
              + mechanical().RhieChowCorrection
                (
                    U(), gradU_, runTime().deltaT()*impKf_
                )
            );

//Info<< sigma() << endl;
            // Under-relaxation the linear system
            //DEqn.relax();

            UEqn.relax();

            // Solve the linear system
            //solverPerfD = DEqn.solve();

            solverPerfU = UEqn.solve();



            // Fixed or adaptive field under-relaxation
            //relaxField(D(), iCorr);
            relaxField(U(), iCorr);



            // Update increment of displacement
            //DD() = D() - D().oldTime();
            DD() = U()*runTime().deltaT();

            D() = D().oldTime() + DD();

            // Update gradient of displacement
            mechanical().grad(D(), gradD());
           // mechanical().grad(U(), gradU());
            gradU_ = fvc::grad(U());

            // Update gradient of displacement increment
            gradDD() = gradD() - gradD().oldTime();
            //gradDD()=gradU()*runTime.deltaT()

            // Calculate the stress using run-time selectable mechanical law
            mechanical().correct(sigma());
//Info<< sigma() << endl;

            // Update impKf to improve convergence
            // Note: impK and rImpK are not updated as they are used for
            // traction boundaries
            //if (iCorr % 10 == 0)
            //{
            //    impKf_ = mechanical().impKf();
            //}
        }
        while
        (
            !converged
            (
                iCorr,
                solverPerfU.initialResidual(), solverPerfU.nIterations(), U()
            )
         && ++iCorr < nCorr()
        );

        // Interpolate cell displacements to vertices
        mechanical().interpolate(D(), pointD());

        // Increment of displacement
        //DD() = D() - D().oldTime();
        //DD()=U()*mag(runTime().deltaT());

        // Increment of point displacement
         pointDD() = pointD() - pointD().oldTime();

        // Velocity
        // U() = fvc::ddt(D());

    }
    while (mesh().update());

    return true;
}

tmp<vectorField> linGeomVelocitySolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch mechanical property
    const scalarField& impK = impK_.boundaryField()[patchID];

    // Patch reciprocal implicit stiffness field
    const scalarField& rImpK = rImpK_.boundaryField()[patchID];

    // Patch gradient
    const tensorField& pGradU = gradU_.boundaryField()[patchID];

    // Patch stress
    const symmTensorField& pSigma = sigma().boundaryField()[patchID];

    // Patch unit normals
    const vectorField n = patch.nf();

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - n*pressure)
                - (n & (pSigma - runTime().deltaT().value()*impK*pGradU))
            )*rImpK/runTime().deltaT().value()
        )
    );
}

}

}

