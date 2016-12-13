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

#include "unsNonLinGeomTotalLagSolid.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "solidTractionFvPatchVectorField.H"
#include "fvcGradf.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(unsNonLinGeomTotalLagSolid, 0);
addToRunTimeSelectionTable
(
    solidModel, unsNonLinGeomTotalLagSolid, dictionary
);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //


scalar unsNonLinGeomTotalLagSolid::residual() const
{
    return
        gMax
        (
            mag(D_.internalField() - D_.prevIter().internalField())
           /max
            (
                gMax(mag(D_.internalField() - D_.oldTime().internalField())),
                SMALL
            )
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

unsNonLinGeomTotalLagSolid::unsNonLinGeomTotalLagSolid(dynamicFvMesh& mesh)
:
    solidModel(typeName, mesh),
    D_
    (
        IOobject
        (
            "D",
            runTime().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    U_
    (
        IOobject
        (
            "U",
            runTime().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimLength/dimTime, vector::zero)
    ),
    pMesh_(mesh),
    pointD_
    (
        IOobject
        (
            "pointD",
            runTime().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pMesh_,
        dimensionedVector("0", dimLength, vector::zero)
    ),
    sigma_
    (
        IOobject
        (
            "sigmaCauchy",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    ),
    sigmaf_
    (
        IOobject
        (
            "sigmaCauchyf",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    ),
    gradD_
    (
        IOobject
        (
            "grad(" + D_.name() + ")",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("0", dimless, tensor::zero)
    ),
    gradDf_
    (
        IOobject
        (
            "grad(" + D_.name() + ")f",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("0", dimless, tensor::zero)
    ),
    F_
    (
        IOobject
        (
            "F",
            runTime().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedTensor("I", dimless, I)
    ),
    Ff_
    (
        IOobject
        (
            "Ff",
            runTime().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedTensor("I", dimless, I)
    ),
    Finv_
    (
        IOobject
        (
            "Finv",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        inv(F_)
    ),
    Finvf_
    (
        IOobject
        (
            "Finvf",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        inv(Ff_)
    ),
    J_
    (
        IOobject
        (
            "J",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        det(F_)
    ),
    Jf_
    (
        IOobject
        (
            "Jf",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        det(Ff_)
    ),
    rho_(mechanical().rho()),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_),
    solutionTol_
    (
        solidProperties().lookupOrDefault<scalar>("solutionTolerance", 1e-06)
    ),
    relativeTol_
    (
        solidProperties().lookupOrDefault<scalar>("relativeTolerance", 1e-03)
    ),
    materialTol_
    (
        solidProperties().lookupOrDefault<scalar>("materialTolerance", 1e-05)
    ),
    infoFrequency_
    (
        solidProperties().lookupOrDefault<int>("infoFrequency", 100)
    ),
    nCorr_(solidProperties().lookupOrDefault<int>("nCorrectors", 10000)),
    g_
    (
        IOobject
        (
            "g",
            runTime().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    nonLinear_(solidProperties().lookupOrDefault<Switch>("nonLinear", true)),
    enforceLinear_(false),
    debug_(solidProperties().lookupOrDefault<Switch>("debug", false)),
    K_
    (
        solidProperties().lookupOrDefault<dimensionedScalar>
        (
            "K",
            dimensionedScalar("K", dimless/dimTime, 0)
        )
    ),
    maxIterReached_(0)
{
    D_.oldTime().oldTime();
    pointD_.oldTime();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> unsNonLinGeomTotalLagSolid::patchDisplacementIncrement
(
    const label patchID
) const
{
    tmp<vectorField> tPatchDisplacement
    (
        new vectorField
        (
            mesh().boundary()[patchID].size(),
            vector::zero
        )
    );
    vectorField& patchDisplacement = tPatchDisplacement();

    patchDisplacement =
        D_.boundaryField()[patchID] - D_.oldTime().boundaryField()[patchID];

    return tPatchDisplacement;
}

void unsNonLinGeomTotalLagSolid::setTraction
(
    const label patchID,
    const vectorField& traction
)
{
    if
    (
        D_.boundaryField()[patchID].type()
     != solidTractionFvPatchVectorField::typeName
    )
    {
        FatalErrorIn("void unsNonLinGeomTotalLagSolid::setTraction(...)")
            << "Boundary condition on " << D_.name()
            <<  " is "
            << D_.boundaryField()[patchID].type()
            << " for patch" << mesh().boundary()[patchID].name()
            << ", instead "
            << solidTractionFvPatchVectorField::typeName
            << abort(FatalError);
    }

    solidTractionFvPatchVectorField& patchD =
        refCast<solidTractionFvPatchVectorField>
        (
            D_.boundaryField()[patchID]
        );

    patchD.traction() = traction;
}


void unsNonLinGeomTotalLagSolid::setPressure
(
    const label patchID,
    const scalarField& pressure
)
{
    if
    (
        D_.boundaryField()[patchID].type()
     != solidTractionFvPatchVectorField::typeName
    )
    {
        FatalErrorIn("void unsNonLinGeomTotalLagSolid::setTraction(...)")
            << "Boundary condition on " << D_.name()
            <<  " is "
            << D_.boundaryField()[patchID].type()
            << " for patch" << mesh().boundary()[patchID].name()
            << ", instead "
            << solidTractionFvPatchVectorField::typeName
            << abort(FatalError);
    }

    solidTractionFvPatchVectorField& patchD =
        refCast<solidTractionFvPatchVectorField>
        (
            D_.boundaryField()[patchID]
        );

    patchD.pressure() = pressure;
}


bool unsNonLinGeomTotalLagSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    int iCorr = 0;
    scalar initialResidual = 0;
    solverPerformance solverPerf;
    scalar res = 1.0;
    scalar maxRes = 0;
    scalar curConvergenceTolerance = solutionTol_;
    solverPerformance::debug = 0;

    do
    {
        if (lduMatrix::debug)
        {
            Info<< "Time: " << runTime().timeName()
                << ", outer iteration: " << iCorr << endl;
        }

        // Store previous iteration to allow under-relaxation and residual
        // calculation
        D_.storePrevIter();

        // Construct momentum equation in total Lagrangian form where gradients
        // are calculated directly at the faces
        fvVectorMatrix DEqn
        (
            rho_*fvm::d2dt2(D_)
         == fvm::laplacian(impKf_, D_, "laplacian(DD,D)")
          - fvc::laplacian(impKf_, D_, "laplacian(DD,D)")
          + fvc::div((Jf_*Finvf_.T() & mesh().Sf()) & sigmaf_)
          + rho_*g_
        );

        // Add damping
        if (K_.value() > SMALL)
        {
            DEqn += K_*rho_*fvm::ddt(D_);
        }

        // Enforce linear to improve convergence
        if (enforceLinear_)
        {
            // Replace nonlinear terms with linear
            // Note: the mechanical law could still be nonlinear
            DEqn +=
                fvc::div((Jf_*Finvf_.T() & mesh().Sf()) & sigmaf_)
              - fvc::div(mesh().Sf() & sigmaf_);
        }

        // Under-relax the linear system
        DEqn.relax();

        // Solve the system
        solverPerf = DEqn.solve();

        // Under-relax displacement field
        D_.relax();

        if (iCorr == 0)
        {
            initialResidual = solverPerf.initialResidual();
        }

        // Interpolate D to pointD
        mechanical().interpolate(D_, pointD_, false);

        // Update gradient of displacement
        mechanical().grad(D_, pointD_, gradD_);
        mechanical().grad(D_, pointD_, gradDf_);

        // Total deformation gradient
        Ff_ = I + gradDf_.T();

        // Inverse of the deformation gradient
        Finvf_ = inv(Ff_);

        // Jacobian of the deformation gradient
        Jf_ = det(Ff_);

        // Calculate the stress using run-time selectable mechanical law
        mechanical().correct(sigmaf_);

        if (nonLinear_ && !enforceLinear_)
        {
            //surfaceScalarField Det = det(I + gradDf_);

            scalar minJf = min(Jf_).value();
            reduce(minJf, minOp<scalar>());

            scalar maxJf = max(Jf_).value();
            reduce(maxJf, maxOp<scalar>());

            if ((minJf < 0.01) || (maxJf > 100))
            {
                Info<< minJf << ", " << maxJf << endl;

                // Enable enforce linear to try improve convergence
                enforceLinear_ = true;
            }
        }

        // Calculate relative momentum residual
        res = residual();

        if (res > maxRes)
        {
            maxRes = res;
        }

        curConvergenceTolerance = maxRes*relativeTol_;

        if (curConvergenceTolerance < solutionTol_)
        {
            curConvergenceTolerance = solutionTol_;
        }

        if (lduMatrix::debug)
        {
            Info<< "Relative residual = " << res << endl;
        }

        if (maxIterReached_ == nCorr_)
        {
            maxIterReached_++;
        }
    }
    while (res > curConvergenceTolerance && ++iCorr < nCorr_);

    // Velocity
    U_ = fvc::ddt(D_);

    // Total deformation gradient
    F_ = I + gradD_.T();

    // Inverse of the deformation gradient
    Finv_ = inv(F_);

    // Jacobian of the deformation gradient
    J_ = det(F_);

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma_);

    // Print summary of residuals
    Info<< solverPerf.solverName() << ": Solving for " << D_.name()
        << ", Initial residual = " << initialResidual
        << ", Final residual = " << solverPerf.initialResidual()
        << ", No outer iterations = " << iCorr << nl
        << " Max relative residual = " << maxRes
        << ", Relative residual = " << res
        << ", enforceLinear = " << enforceLinear_ << endl;

    lduMatrix::debug = 1;

    if (nonLinear_ && enforceLinear_)
    {
        return false;
    }

    return true;
}


tmp<vectorField> unsNonLinGeomTotalLagSolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch mechanical property
    const scalarField& impK = impKf_.boundaryField()[patchID];

    // Patch reciprocal implicit stiffness field
    const scalarField& rImpK = rImpK_.boundaryField()[patchID];

    // Patch gradient
    const tensorField& gradD = gradDf_.boundaryField()[patchID];

    // Patch stress
    const symmTensorField& sigma = sigmaf_.boundaryField()[patchID];

    // Patch total deformation gradient inverse
    const tensorField& Finv = Finvf_.boundaryField()[patchID];

    // Patch total Jacobian
    const scalarField& J = Jf_.boundaryField()[patchID];

    // Patch unit normals (initial configuration)
    const vectorField n = patch.nf();

    // Patch unit normals (deformed configuration)
    const vectorField nCurrent = J*Finv.T() & n;

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - n*pressure)
              - (nCurrent & sigma)
              + (n & (impK*gradD))
            )*rImpK
        )
    );
}


void unsNonLinGeomTotalLagSolid::updateTotalFields()
{
    mechanical().updateTotalFields();
}


void unsNonLinGeomTotalLagSolid::writeFields(const Time& runTime)
{
    // Calculate equivalent (von Mises) stress
    volScalarField sigmaEq
    (
        IOobject
        (
            "sigmaCauchyEq",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt((3.0/2.0)*magSqr(dev(sigma_)))
    );

    Info<< "Max sigmaCauchyEq = " << gMax(sigmaEq) << endl;

    solidModel::writeFields(runTime);
}


void unsNonLinGeomTotalLagSolid::end()
{
    if (maxIterReached_ > 0)
    {
        WarningIn(type() + "::end()")
            << "The maximum momentum correctors were reached in "
            << maxIterReached_ << " time-steps" << nl << endl;
    }
    else
    {
        Info<< "The momentum equation converged in all time-steps"
            << nl << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
