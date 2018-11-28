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

#include "linGeomTotalDispSolidDamage.H"
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

defineTypeNameAndDebug(linGeomTotalDispSolidDamage, 0);
addToRunTimeSelectionTable(physicsModel, linGeomTotalDispSolidDamage, solid);
addToRunTimeSelectionTable(solidModel, linGeomTotalDispSolidDamage, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linGeomTotalDispSolidDamage::linGeomTotalDispSolidDamage
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_),
    damage_
    (
        IOobject
        (
            "damage",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    )
        //    mu
//    (
//        IOobject
//        (
//            "interpolate(mu)",
//            runTime.timeName(),
//            mesh(),
//            IOobject::NO_READ,
//            IOobject::NO_WRITE
//        ),
//        mesh(),
//        dimensionedScalar("0", dimPressure, 0.0)
//    ),
//    lambda
//    (
//        IOobject
//        (
//            "interpolate(lambda)",
//            runTime.timeName(),
//            mesh(),
//            IOobject::NO_READ,
//            IOobject::NO_WRITE
//        ),
//        mesh(),
//        dimensionedScalar("0", dimPressure, 0.0)
//    )
//    mu
//    (
//        IOobject
//        (
//            "mu",
//            runTime.timeName(),
//            mesh(),
//            IOobject::NO_READ,
//            IOobject::NO_WRITE
//        ),
//        mesh(),
//        dimensionedScalar("0", dimPressure, 0.0)
//    ),
//    lambda
//    (
//        IOobject
//        (
//            "lambda",
//            runTime.timeName(),
//            mesh(),
//            IOobject::NO_READ,
//            IOobject::NO_WRITE
//        ),
//        mesh(),
//        dimensionedScalar("0", dimPressure, 0.0)
//    )
//solidModelDict().lookupOrDefault<scalar>("solutionTolerance", 1e-06)
//    S0_(solidModelDict.lookup("S0")),#
//    S0_(dict.lookup("sD")),
//    h_(dict.lookup("h")),
//    b_(dict.lookup("b")),
//    epsilonD_(dict.lookup("epsilonD"))

{
//const PtrList<mechanicalLaw>& mechLaws = mechanical();
//    const linearElasticMisesPlastic& mech = refCast<const linearElasticMisesPlastic>(mechLaws[0]);
//    mu = mech.mu();
//    lambda = mech.lambda();
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool linGeomTotalDispSolidDamage::evolve()
{
    Info<< "Evolving solid solver" << endl;

    // Mesh update loop
    do
    {
        int iCorr = 0;
        lduSolverPerformance solverPerfD;
        blockLduMatrix::debug = 0;

        Info<< "Solving the momentum equation for D" << endl;

        // Stablisation viscosity
        // const dimensionedScalar eta_ =
        //    solidModelDict().lookupOrDefault<dimensionedScalar>
        //    (
        //        "numericalViscosity",
        //        dimensionedScalar("eta", dimless/dimTime, 0.0)
        //    );

        // Momentum equation loop
        do
        {
            // Store fields for under-relaxation and residual calculation
            D().storePrevIter();

            // Linear momentum equation total displacement form
            fvVectorMatrix DEqn
            (
                rho()*fvm::d2dt2(D())
                //+ eta_*rho_*fvm::ddt(D())
             == fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
              - fvc::laplacian(impKf_, D(), "laplacian(DD,D)")
              + fvc::div(sigma(), "div(sigma)")
              + rho()*g()
              + mechanical().RhieChowCorrection(D(), gradD())
            );

            // Under-relaxation the linear system
            DEqn.relax();

            // Solve the linear system
            solverPerfD = DEqn.solve();

            // Fixed or adaptive field under-relaxation
            relaxField(D(), iCorr);

            // Update increment of displacement
            DD() = D() - D().oldTime();

            // Update gradient of displacement
            mechanical().grad(D(), gradD());

            // Update gradient of displacement increment
            gradDD() = gradD() - gradD().oldTime();

            // Calculate the stress using run-time selectable mechanical law
            mechanical().correct(sigma());

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
                solverPerfD.initialResidual(), solverPerfD.nIterations(), D()
            )
         && ++iCorr < nCorr()
        );

        // Interpolate cell displacements to vertices
        mechanical().interpolate(D(), pointD());

        // Increment of displacement
        DD() = D() - D().oldTime();

        // Increment of point displacement
        pointDD() = pointD() - pointD().oldTime();

        // Velocity
        U() = fvc::ddt(D());
    }
    while (mesh().update());


    // Calculate Lemaitre damage field

    dimensionedScalar S0_=dimensionedScalar(solidProperties().lookup("s0"));
    dimensionedScalar h_=dimensionedScalar(solidProperties().lookup("h"));
    dimensionedScalar b_=readScalar(solidProperties().lookup("b"));

    scalar mu = readScalar(solidProperties().lookup("mu"));
    scalar lambda = readScalar(solidProperties().lookup("lambda"));
    dimensionedScalar epsilonD_ =
        dimensionedScalar(solidProperties().lookup("epsilonD"));

    // Lookup Kirchhoff stress field
    const volSymmTensorField& sigma =
        mesh().lookupObject<volSymmTensorField>("sigma");

    // Calculate equivalent Cauchy stress
    const volScalarField sigmaEq =
        sqrt((3.0/2.0)*magSqr(dev(sigma)))
        + dimensionedScalar("SMALL", dimPressure, SMALL);

    // Calculate hydrostatic Cauchy stress
    const volScalarField sigmaHyd = (1.0/3.0)*tr(sigma);

    // Calculate constraint/triality
    const volScalarField triaxiality = sigmaHyd/sigmaEq;

    // Lookup plastic strain increment field
    const volScalarField& DEpsilonPEq =
        mesh().lookupObject<volScalarField>("DEpsilonPEq");

    // Lookup plastic strain field
    const volScalarField& epsilonPEq =
        mesh().lookupObject<volScalarField>("epsilonPEq");

    // Young's modulus
    const scalar E = mu*(3.0*lambda + 2.0*mu)/(lambda + mu);

    // Poisson's ration
    const  scalar nu = lambda/(2*(lambda + mu));
    volScalarField Y
        (
            IOobject
            (
                "Y",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimPressure, 0.0),
            zeroGradientFvPatchScalarField::typeName
        );

    volScalarField Ystar
        (
            IOobject
            (
                "Ystar",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimless, 0.0),
            zeroGradientFvPatchScalarField::typeName
        );

    scalarField& YI = Y.internalField();
    scalarField& YstarI = Ystar.internalField();
    const scalarField& sigmaEqI = sigmaEq.internalField();
    const scalarField& epsilonPEqI = epsilonPEq.internalField();
    const scalarField& triaxialityI = triaxiality.internalField();
    const scalar epsilonDValue = epsilonD_.value();
    scalarField& damageI = damage_.internalField();
    const scalarField damageOldI = damage_.oldTime().internalField();
    const scalarField& DepsilonPeqI = DEpsilonPEq.internalField();

    if (debug)
    {
        Info<< nl << "Calculating damage" << endl;
    }

    forAll(YI, cellI)
    {
        if (debug > 1)
        {
            Info<< nl << "Cell = " << cellI << endl;
        }

        if
        (
            triaxialityI[cellI] > (-1.0/3.0)
         && epsilonPEqI[cellI] > epsilonDValue
        )
        {
            const scalar denom =
                (2.0/3.0)*(1.0 + nu)
                + 3.0*(1.0 - 2.0*nu)*pow(triaxialityI[cellI], 2.0);

            YI[cellI]=
                -(
                    pow(sigmaEqI[cellI], 2.0)
                    /(
                        2.0*E
                        // *pow(1.0 - damageLemaitreI[cellI], 2.0)
                    )
                )*denom;

            if (debug > 1)
            {
                Info <<"YI: " << YI[cellI] << endl;
            }

            YstarI[cellI] =
                pow(-YI[cellI]/S0_.value(), b_.value())*DepsilonPeqI[cellI];

            if (debug > 1)
            {
                Info<< "YstarI: " << YstarI[cellI] << endl;
            }

            scalar d1 = damageI[cellI];

            for (int i = 1; i <= 100; ++i)
            {
                scalar fx = d1 - YstarI[cellI]*pow(1.0-d1, -2.0*b_.value());

                if (debug > 1)
                {
                    Info<< "Newton iteration = " << i << nl
                        << "fx: " << fx << endl;
                }
                scalar dfx =
                    1
                  + YstarI[cellI]*2.0*b_.value()*pow(1.0-d1, -2.0*b_.value()-1);

                if (debug > 1)
                {
                    Info<< "dfx: " << dfx << endl;
                }

                d1 = d1 - fx/dfx;

                if (debug > 1)
                {
                    Info<< "d1: " << d1 << endl;
                }
            }

            damageI[cellI] = d1;
        }
    }

    Y.correctBoundaryConditions();
    Ystar.correctBoundaryConditions();
    damage_.correctBoundaryConditions();

    volScalarField term1("ystar", Ystar);
    term1.write();

    return true;
}


tmp<vectorField> linGeomTotalDispSolidDamage::tractionBoundarySnGrad
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
    const tensorField& pGradD = gradD().boundaryField()[patchID];

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
              - (n & (pSigma - impK*pGradD))
            )*rImpK
        )
    );
}



//        forAll(YI, cellI)
//        {
//            // Calculate damage
//            if
//            (
//                triaxialityI[cellI] > (-1.0/3.0)
//             && epsilonPEqI[cellI] > epsilonDValue
//            )
//            {
//                // Calculate Lemaitre Y parameter
//                const scalar denom =
//                    (2.0/3.0)*(1.0 + nu)
//                  + 3.0*(1.0 - 2.0*nu)*pow(triaxialityI[cellI], 2.0);

//Info <<"denom: "<< denom << endl;

//                if (triaxialityI[cellI] > 0.0)
//                {
//                    // Damage in tension
//                    YI[cellI]=
//                       -(
//                            pow(sigmaEqI[cellI], 2.0)
//                           /(
//                               2.0*E
//                              *pow(1.0 - damageLemaitreI[cellI], 2.0)
//                           )
//                        )*denom;
//Info << "Y: "<<YI[cellI] << endl;
//                }
//                else
//                {
//                    // Damage in compression is less than in tension; this
//                    // depends on the h_ parameter
//                    YI[cellI] =
//                        -(
//                            h_.value()*pow(sigmaEqI[cellI], 2.0)
//                           /(
//                                2.0*E
//                               *pow
//                                (
//                                    1.0 - h_.value()*damageLemaitreI[cellI],
//                                    2.0
//                                )
//                            )
//                        )*denom;
//Info << "compressionY: " <<YI[cellI] << endl;
//
//                }

//                YstarI[cellI]=pow(-YI[cellI]/S0_.value(), b_.value());




//Info << "Ystar: " << YstarI[cellI] <<endl;
//            }
//        }

//     Y.correctBoundaryConditions();
//     Ystar.correctBoundaryConditions();
//     damage_ = damage_.oldTime() + Ystar*DEpsilonPEq/(1 - damage_);
//volScalarField term1("ystar", Ystar);
//term1.write();


//Info <<damage_<<endl;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
