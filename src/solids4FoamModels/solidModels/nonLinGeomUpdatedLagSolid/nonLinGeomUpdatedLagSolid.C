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

#include "nonLinGeomUpdatedLagSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(nonLinGeomUpdatedLagSolid, 0);
addToRunTimeSelectionTable
(
    physicsModel, nonLinGeomUpdatedLagSolid, solid
);
addToRunTimeSelectionTable
(
    solidModel, nonLinGeomUpdatedLagSolid, dictionary
);


// * * * * * * * * * *  Private Member Functions  * * * * * * * * * * * * * * //

void nonLinGeomUpdatedLagSolid::updateTangentModulus
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

nonLinGeomUpdatedLagSolid::nonLinGeomUpdatedLagSolid
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
    relF_
    (
        IOobject
        (
            "relF",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        I + gradDD().T()
    ),
    relFinv_
    (
        IOobject
        (
            "relFinv",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        inv(relF_)
    ),
    relJ_
    (
        IOobject
        (
            "relJ",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        det(relF_)
    ),
    rho_
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mechanical().rho()
    ),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool nonLinGeomUpdatedLagSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    int iCorr = 0;
    lduSolverPerformance solverPerfDD;
    blockLduMatrix::debug = 0;

    Info<< "Solving the updated Lagrangian form of the momentum equation for DD"
        << endl;

    const Switch numericalTangent =
        solidProperties().lookupOrDefault<Switch>("numericalTangent", false);
    Info<< "numericalTangent: " << numericalTangent << endl;

    // Momentum equation loop
    do
    {
        // Store fields for under-relaxation and residual calculation
        DD().storePrevIter();

        // Update implicit stiffness (tangent matrix) to improve convergence
        if (numericalTangent)
        {
            updateTangentModulus
            (
                impK_, gradDD(), relF_, relJ_, relFinv_, sigma()
            );
        }
        else if (mechanical().updateImpK())
        {
            impK_ = mechanical().impK();
        }

        // Update surface and reciprocal fields
        impKf_ = fvc::interpolate(impK_);
        rImpK_ = 1.0/impK_;

        // Momentum equation incremental updated Lagrangian form
        fvVectorMatrix DDEqn
        (
            fvm::d2dt2(rho_, DD())
          + fvc::d2dt2(rho_.oldTime(), D().oldTime())
         == fvm::laplacian(impKf_, DD(), "laplacian(DDD,DD)")
          - fvc::laplacian(impKf_, DD(), "laplacian(DDD,DD)")
          + fvc::div(relJ_*relFinv_ & sigma(), "div(sigma)")
          + rho_*g()
          //+ mechanical().RhieChowCorrection(DD(), gradDD())
          + mechanical().RhieChowCorrection(DD(), gradDD(), impKf_)
        );

        // Under-relax the linear system
        DDEqn.relax();

        // Solve the linear system
        solverPerfDD = DDEqn.solve();

        // Under-relax the DD field using fixed or adaptive under-relaxation
        relaxField(DD(), iCorr);

        // Update the total displacement
        D() = D().oldTime() + DD();

        // Update gradient of displacement increment
        mechanical().grad(DD(), gradDD());

        // Relative deformation gradient
        relF_ = I + gradDD().T();

        // Inverse relative deformation gradient
        relFinv_ = inv(relF_);

        // Total deformation gradient
        F_ = relF_ & F_.oldTime();

        // Relative Jacobian (Jacobian of relative deformation gradient)
        relJ_ = det(relF_);

        // Jacobian of deformation gradient
        J_ = relJ_*J_.oldTime();

        // Update the momentum equation inverse diagonal field
        // This may be used by the mechanical law when calculating the
        // hydrostatic pressure
        const volScalarField DEqnA("DEqnA", DDEqn.A());

        // Calculate the stress using run-time selectable mechanical law
        mechanical().correct(sigma());
    }
    while
    (
       !converged
        (
            iCorr,
            solverPerfDD.initialResidual(),
            solverPerfDD.nIterations(),
            DD()
        )
     && ++iCorr < nCorr()
    );

    // Update gradient of total displacement
    gradD() = fvc::grad(D().oldTime() + DD());

    // Total displacement
    D() = D().oldTime() + DD();

    // Update pointDD as it used by FSI procedure
    mechanical().interpolate(DD(), pointDD());

    // Total displacement at points
    pointD() = pointD().oldTime() + pointDD();

    // Velocity
    U() = fvc::ddt(D());

    return true;
}


tmp<vectorField> nonLinGeomUpdatedLagSolid::tractionBoundarySnGrad
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

    // Patch gradient
    const tensorField& pGradDD = gradDD().boundaryField()[patchID];

    // Patch Cauchy stress
    const symmTensorField& pSigma = sigma().boundaryField()[patchID];

    // Patch relative deformation gradient inverse
    const tensorField& relFinv = relFinv_.boundaryField()[patchID];

    // Patch unit normals (updated configuration)
    const vectorField n = patch.nf();

    // Patch unit normals (deformed configuration)
    vectorField nCurrent = relFinv.T() & n;
    nCurrent /= mag(nCurrent);

    // Testing: let us instead calculate the deformed normals by interpolating
    // displacements to the points and calculating the normals on the deformed
    // patch; as this is how we will actually move the mesh, it will be more
    // consistent.
    // This, however, begs the question: is the cell-centred deformation
    // gradient field 'F' consistent with our point displacement field?"
    // i.e. we can calculate the deformed cell volumes two ways (at least):
    //     1. V = J*Vold
    //     2. Move the mesh with pointD and then directly calculate V
    // The answers from 1. and 2. are only approximately equal: this causes a
    // slight inconsistency. The equalavent can be said for the deformed face
    // areas.
    // In Maneeratana, the mesh is never moved, instead method 1. is used for
    // the deformed volumes and areas.

    // standAlonePatch deformedPatch =
    //     standAlonePatch
    //     (
    //         mesh().boundaryMesh()[patchID].localFaces(),
    //         mesh().boundaryMesh()[patchID].localPoints()
    //     );

    // // Calculate the deformed points
    // const pointField deformedPoints =
    //     mechanical().volToPoint().interpolate
    //     (
    //         mesh().boundaryMesh()[patchID],
    //         DD_
    //     )
    //   + mesh().boundaryMesh()[patchID].localPoints();

    // // Move the standAlonePatch points
    // const_cast<pointField&>(deformedPatch.points()) = deformedPoints;

    // // Patch unit normals (deformed configuration)
    // const vectorField& nCurrent = deformedPatch.faceNormals();

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - nCurrent*pressure)
              - (nCurrent & pSigma)
              + impK*(n & pGradDD)
            )*rImpK
        )
    );
}


void nonLinGeomUpdatedLagSolid::updateTotalFields()
{
    // Density
    rho_ = rho_.oldTime()/relJ_;

    // Move the mesh to the deformed configuration
    const vectorField oldPoints = mesh().allPoints();
    moveMesh(oldPoints, DD(), pointDD());

    solidModel::updateTotalFields();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
