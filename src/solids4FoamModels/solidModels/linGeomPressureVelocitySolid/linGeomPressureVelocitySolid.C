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

#include "linGeomPressureVelocitySolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "linearElastic.H"
#include "findRefCell.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(linGeomPressureVelocitySolid, 0);
addToRunTimeSelectionTable(physicsModel, linGeomPressureVelocitySolid, solid);
addToRunTimeSelectionTable(solidModel, linGeomPressureVelocitySolid, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linGeomPressureVelocitySolid::linGeomPressureVelocitySolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_),
    gradU_(fvc::grad(U())),
    p_
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
    gradp_(fvc::grad(p_)),
    mu_
    (
        IOobject
        (
            "mu",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("0", dimPressure, 0.0)
    ),
    muf_
    (
        IOobject
        (
            "interpolate(mu)",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("0", dimPressure, 0.0)
    ),
    lambdaf_
    (
        IOobject
        (
            "interpolate(lambda)",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("0", dimPressure, 0.0)
    ),
    kf_
    (
        IOobject
        (
            "k",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("0", dimPressure, 0.0)
    ),
currentRho_(rho())


{    // Cast the mechanical law to a linearElastic mechanicalLaw
const PtrList<mechanicalLaw>& mechLaws = mechanical();
    const linearElastic& mech = refCast<const linearElastic>(mechLaws[0]);

    // Set mu and lambda fields
    mu_ = mech.mu();
    muf_ = mech.mu();
    //muf_ = fvc::interpolate(mu_);
    lambdaf_ = mech.lambda();
    kf_=mech.K();
    p_.oldTime();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool linGeomPressureVelocitySolid::evolve()
{
    //nICorr_(solidProperties().lookupOrDefault<int>("nInnerCorrectors", 3))
    const int nICorr_ = readInt(solidProperties().lookup("nInnerCorrectors"));
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p_, solidProperties(), pRefCell, pRefValue);

    //Momentum

    //volScalarField& rho1 = rho();
    //const volScalarField& rho1 = rho();


    //const int nICorr =
    //    readInt(solidProperties().lookupOrDefault("nInnerCorrectors", 3));



    for (int oCorr = 0; oCorr < nCorr(); oCorr++)
    {
        Info<< "oCorr: " << oCorr << endl;

        U().storePrevIter();

        fvVectorMatrix UEqn
        (
            fvm::ddt(currentRho_, U())
            //fvm::ddt(rho(), U())
          - fvm::laplacian(2*muf_*runTime().deltaT(), U())
            //+fvc::div(2*muf_*runTime().deltaT(), fvc::grad(U()))
          + fvc::div(2*mu_*runTime().deltaT()*gradU_)
          - fvc::div(dev(sigma()))
        );

        UEqn.relax();

        solve(UEqn == -gradp_);

        //UEqn.solve();
        //Info << U() << endl;

        volScalarField rAU = 1/UEqn.A();


        // Start of Piso Loop
        for (int iCorr = 0; iCorr < nICorr_; iCorr++)
        {
            p_.storePrevIter();

            fvScalarMatrix PEqn
            (
                fvm::ddt(rho()/kf_, p_)
                //+fvc::div(rho(), UEqn.H())
              + fvc::div(rho()*rAU*UEqn.H())
                //+ fvc::div(((rho()*rAU)/kf_)*UEqn.H()*p_)
              + fvm::div
                (
                    mesh().Sf()
                  & fvc::interpolate(((rho()*rAU)/kf_)*UEqn.H()),
                    p_
                )
                //-fvm::laplacian(currentRho_*rAU, p_)
             == fvm::laplacian(currentRho_*rAU, p_)
                // Rhie-Chow
              + (
                    fvc::laplacian(currentRho_*rAU, p_)
                  - fvc::div(currentRho_*rAU*gradp_)
                )
            );

            //           fvScalarMatrix PEqn
            //           (
            //           fvm::ddt(rho()/kf_, p_)
            //           +fvm::div((rho()/kf_), U()*p_)
            //           ==
            //           -fvc::div(rho()*U())
            //           );

            PEqn.setReference(pRefCell, pRefValue);

            PEqn.relax();

            PEqn.solve();
            //solve(PEqn==0);

            p_.relax();

            gradp_ = fvc::grad(p_);

            U() = UEqn.H()*rAU - rAU*gradp_;
            U().correctBoundaryConditions();

            gradU_ = fvc::grad(U());

            // currentRho_ = rho()*(1 + (p_/kf_));
            currentRho_ = rho();
            //Info << currentRho_ << endl;

            //phi() = (fvc::interpolate(U()) & mesh.Sf());
        }

        DD() = U()*runTime().deltaT();

        D()= DD() + D().oldTime();

        mechanical().grad(D(), gradD());
        // mechanical().grad(U(), gradU());

        // Update gradient of displacement increment
        gradDD() = gradD() - gradD().oldTime();
        //gradDD()=gradU()*runTime.deltaT()

        // Calculate the stress using run-time selectable mechanical law
        mechanical().correct(sigma());

        // Replace the volumetric component of stress
        sigma() = dev(sigma()) - p_*I;
    }

    return true;
}

tmp<vectorField> linGeomPressureVelocitySolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch mechanical property
    const scalarField impK = 2.0*muf_.boundaryField()[patchID];

    // Patch reciprocal implicit stiffness field
    const scalarField rImpK = 1.0/impK;

    // Patch gradient
    const tensorField& pGradU = gradU_.boundaryField()[patchID];

    // Patch stress
    const symmTensorField& pSigma = sigma().boundaryField()[patchID];

    // Patch pressure
    const scalarField& pP = p_.boundaryField()[patchID];

    // Patch unit normals
    const vectorField n = patch.nf();

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - n*pressure)
              - (n & (dev(pSigma) - runTime().deltaT().value()*impK*pGradU))
              + n*pP
            )*rImpK/runTime().deltaT().value()
        )
    );
}

}

}
