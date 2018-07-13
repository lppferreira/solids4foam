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
            "interpolate(lambda)",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("0", dimPressure, 0.0)
    )


{    // Cast the mechanical law to a linearElastic mechanicalLaw
const PtrList<mechanicalLaw>& mechLaws = mechanical();
    const linearElastic& mech = refCast<const linearElastic>(mechLaws[0]);

    // Set mu and lambda fields
    muf_ = mech.mu();
    lambdaf_ = mech.lambda(); 
    kf_=mech.K();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool linGeomPressureVelocitySolid::evolve()
{
 
//Momentum

volScalarField rho1=rho();

const int nICorr = readInt(solidProperties().lookupOrDefault("nInnerCorrectors", 3));


for (int oCorr = 0; oCorr < nCorr(); oCorr++)
{


fvVectorMatrix UEqn

(
fvm::ddt(rho1, U())
-fvm::laplacian(2*muf_*runTime().deltaT(), U())
+fvc::div(2*muf_*runTime().deltaT(), fvc::grad(U()))
-fvc::div(dev(sigma))
);

solve(UEqn==-fvc::grad(p_));

volScalarField rAU=1/UEqn.A();

       for (int iCorr = 0; iCorr < nICorr; iCorr++)
       //Start of Piso Loop
       {

       fvVectorMatrix PEqn
       (
       fvm::ddt(rho()/k, p_)
       +fvc::div(rho(), UEqn.H()*rAU)
       +fvc::div(rho()/k, UEqn.H()*rAU)
       -fvm::laplacian(rho1*rAU, p_)
       );

       solve(PEqn==0);

       U() = UEqn.H()*rAU-rAU*fvm::grad(p_);

       rho1=rho()*(1+(p_/k));

       //phi() = (fvc::interpolate(U()) & mesh.Sf());

       }

DD()=U()*runTime().deltaT();

D()=D().oldTime+U()*runTime().deltaT();

mechanical().grad(D(), gradD());
// mechanical().grad(U(), gradU());

// Update gradient of displacement increment
gradDD() = gradD() - gradD().oldTime();  //gradDD()=gradU()*runTime.deltaT()

// Calculate the stress using run-time selectable mechanical law
mechanical().correct(sigma());

gradU_ = fvc::grad(U());

}

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

