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

#include "nonLinGeomTotalLagTotalDispSolid.H"
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

defineTypeNameAndDebug(nonLinGeomTotalLagTotalDispSolid, 0);
addToRunTimeSelectionTable
(
    physicsModel, nonLinGeomTotalLagTotalDispSolid, solid
);
addToRunTimeSelectionTable
(
    solidModel, nonLinGeomTotalLagTotalDispSolid, dictionary
);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //


void nonLinGeomTotalLagTotalDispSolid::predict()
{
    Info << "Predicting solid model" << endl;

    // Predict D using the velocity field
    // Note: the case may be steady-state but U can still be calculated using a
    // transient method
    D() = D().oldTime() + U()*runTime().deltaT();

    // Update gradient of displacement
    mechanical().grad(D(), gradD());

    // Total deformation gradient
    F_ = I + gradD().T();

    // Inverse of the deformation gradient
    Finv_ = inv(F_);

    // Jacobian of the deformation gradient
    J_ = det(F_);

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());
}


void nonLinGeomTotalLagTotalDispSolid::updateTangentModulus
(
    volScalarField& impK,
    volTensorField& gradD,
    volTensorField& F,
    volScalarField& J,
    volTensorField& Finv,
    volSymmTensorField& sigma
)
{
    // Take a copy of F and sigma as we will reset their values at the end
    const volTensorField unperturbedGradD("unperturbedGradD", gradD);
    const volTensorField unperturbedF("unperturbedFinv", F);
    const volScalarField unperturbedJ("unperturbedJ", J);
    const volTensorField unperturbedFinv("unperturbedFinv", Finv);
    const volSymmTensorField unperturbedSigma("unperturbedSigma", sigma);

    // Use numerical differentiation to approx the diagonal components of the
    // tangent matrix
    // C1111 = (sigma(F + deltaF11) - sigma(F))/deltaF11;
    // C2222 = (sigma(F + deltaF22) - sigma(F))/deltaF22;
    // C3333 = (sigma(F + deltaF33) - sigma(F))/deltaF33;

    // Delta for numerical differentiation
    //const scalar deltaGradD = sqrt(SMALL);
    scalar deltaGradD = sqrt(SMALL);
    //const scalar deltaGradD = 0.1;

    // Mechanical laws look gradD or gradDD

    // Calculate C1111 by perturbing F11/gradD11
    gradD = unperturbedGradD + tensor(deltaGradD, 0, 0, 0, 0, 0, 0, 0, 0);
    F = I + gradD.T();
    J = det(F);
    Finv = inv(F);
    mechanical().correct(sigma);
    const volScalarField C1111 =
        (
            (
                (J*Finv & sigma)
              - (unperturbedJ*unperturbedFinv & unperturbedSigma)
            )/deltaGradD
        )().component(tensor::XX);

    // Calculate C2222 by perturbing F22/gradD22
    gradD = unperturbedGradD + tensor(0, 0, 0, 0, deltaGradD, 0, 0, 0, 0);
    F = I + gradD.T();
    J = det(F);
    Finv = inv(F);
    mechanical().correct(sigma);
    const volScalarField C2222 =
        (
            (
                (J*Finv & sigma)
              - (unperturbedJ*unperturbedFinv & unperturbedSigma)
            )/deltaGradD
        )().component(tensor::YY);

    // Calculate C3333 by perturbing F33/gradD33
    gradD = unperturbedGradD + tensor(0, 0, 0, 0, 0, 0, 0, 0, deltaGradD);
    F = I + gradD.T();
    J = det(F);
    Finv = inv(F);
    mechanical().correct(sigma);
    const volScalarField C3333 =
        (
            (
                (J*Finv & sigma)
              - (unperturbedJ*unperturbedFinv & unperturbedSigma)
            )/deltaGradD
        )().component(tensor::ZZ);

    // Set impK as the average (or max) of the diagonal components
    impK = (C1111 + C2222 + C3333)/3.0;
    //impK = max(C1111, max(C2222, C3333));
    //impK = min(C1111, min(C2222, C3333));

    // Reset gradD and sigma
    gradD = 1.0*unperturbedGradD;
    F = 1.0*unperturbedF;
    J = 1.0*unperturbedJ;
    Finv = 1.0*unperturbedFinv;
    sigma = 1.0*unperturbedSigma;

    // Reset stress in the mechanical law
    mechanical().correct(sigma);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nonLinGeomTotalLagTotalDispSolid::nonLinGeomTotalLagTotalDispSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    F_
    (
        IOobject
        (
            "F",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedTensor("I", dimless, I)
    ),
    Finv_
    (
        IOobject
        (
            "Finv",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        inv(F_)
    ),
    J_
    (
        IOobject
        (
            "J",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        det(F_)
    ),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_),
    predict_(solidProperties().lookupOrDefault<Switch>("predict", false))
{
    Info<< type() << nl
        << "    predict: " << predict_ << nl
        << "    updateImpK: " << Switch(mechanical().updateImpK()) << endl;

    Info<< "initial impK:" << nl
        << "    min: " << min(impK_.internalField()) << nl
        << "    max: " << min(impK_.internalField()) << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool nonLinGeomTotalLagTotalDispSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    if (predict_)
    {
        predict();
    }

    int iCorr = 0;
    lduSolverPerformance solverPerfD;
    blockLduMatrix::debug = 0;

    Info<< "Solving the total Lagrangian form of the momentum equation for D"
        << endl;

    const Switch numericalTangent =
        solidProperties().lookupOrDefault<Switch>("numericalTangent", false);
    Info<< "numericalTangent: " << numericalTangent << endl;

    // Momentum equation loop
    do
    {
        // Store fields for under-relaxation and residual calculation
        D().storePrevIter();

        // Update implicit stiffness (tangent matrix) to improve convergence
        if (numericalTangent)
        {
            updateTangentModulus(impK_, gradD(), F_, J_, Finv_, sigma());
        }
        else if (mechanical().updateImpK())
        {
            impK_ = mechanical().impK();
        }

        // Update surface and reciprocal fields
        impKf_ = fvc::interpolate(impK_);
        rImpK_ = 1.0/impK_;

        // Momentum equation total displacement total Lagrangian form
        fvVectorMatrix DEqn
        (
            rho()*fvm::d2dt2(D())
         == fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
          - fvc::laplacian(impKf_, D(), "laplacian(DD,D)")
          + fvc::div(J_*Finv_ & sigma(), "div(sigma)")
          + rho()*g()
          + mechanical().RhieChowCorrection(D(), gradD(), impKf_)
        );

        // Under-relax the linear system
        DEqn.relax();

        // Solve the linear system
        solverPerfD = DEqn.solve();

        // Fixed or adaptive field under-relaxation
        relaxField(D(), iCorr);

        // Increment of displacement
        DD() = D() - D().oldTime();

        // Update gradient of displacement
        mechanical().grad(D(), gradD());

        // Update gradient of displacement increment
        gradDD() = gradD() - gradD().oldTime();

        // Total deformation gradient
        F_ = I + gradD().T();

        // Inverse of the deformation gradient
        Finv_ = inv(F_);

        // Jacobian of the deformation gradient
        J_ = det(F_);

        // Update the momentum equation inverse diagonal field
        // This may be used by the mechanical law when calculating the
        // hydrostatic pressure
        const volScalarField DEqnA("DEqnA", DEqn.A());

        // Calculate the stress using run-time selectable mechanical law
        mechanical().correct(sigma());
    }
    while
    (
       !converged
        (
            iCorr,
            solverPerfD.initialResidual(),
            solverPerfD.nIterations(),
            D()
        ) && ++iCorr < nCorr()
    );

    // Interpolate cell displacements to vertices
    mechanical().interpolate(D(), pointD());

    // Increment of point displacement
    pointDD() = pointD() - pointD().oldTime();

    // Velocity
    U() = fvc::ddt(D());

    return true;
}


tmp<vectorField> nonLinGeomTotalLagTotalDispSolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch implicit stiffness field
    const scalarField& impK = impK_.boundaryField()[patchID];

    // Patch reciprocal implicit stiffness field
    const scalarField& rImpK = rImpK_.boundaryField()[patchID];

    // TEST
    // const scalarField impK = mechanical().impK()().boundaryField()[patchID];
    // const scalarField rImpK = 1.0/impK;

    // Patch gradient
    const tensorField& pGradD = gradD().boundaryField()[patchID];

    // Patch Cauchy stress
    const symmTensorField& pSigma = sigma().boundaryField()[patchID];

    // Patch total deformation gradient inverse
    const tensorField& Finv = Finv_.boundaryField()[patchID];

    // Patch unit normals (initial configuration)
    const vectorField n = patch.nf();

    // Patch unit normals (deformed configuration)
    vectorField nCurrent = Finv.T() & n;
    nCurrent /= mag(nCurrent);

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - nCurrent*pressure)
              - (nCurrent & pSigma)
              + impK*(n & pGradD)
            )*rImpK
        )
    );
}


void nonLinGeomTotalLagTotalDispSolid::writeFields(const Time& runTime)
{
    solidModel::writeFields(runTime);

    // Info<< "Writing impK" << nl
    //     << "    " << min(impK_) << nl
    //     << "    " << max(impK_) << nl
    //     << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
