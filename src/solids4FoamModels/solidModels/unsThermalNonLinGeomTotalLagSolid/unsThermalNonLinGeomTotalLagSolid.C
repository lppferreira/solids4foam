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

#include "unsThermalNonLinGeomTotalLagSolid.H"
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

defineTypeNameAndDebug(unsThermalNonLinGeomTotalLagSolid, 0);
addToRunTimeSelectionTable
(
    physicsModel, unsThermalNonLinGeomTotalLagSolid, solid
);
addToRunTimeSelectionTable
(
    solidModel, unsThermalNonLinGeomTotalLagSolid, dictionary
);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

bool unsThermalNonLinGeomTotalLagSolid::converged
(
    const int iCorr,
    const lduSolverPerformance& solverPerfD,
    const lduSolverPerformance& solverPerfT,
    const volVectorField& D,
    const volScalarField& T
)
{
    // We will check a number of different residuals for convergence
    bool converged = false;

    // Calculate relative residuals
    const scalar absResidualT =
        gMax(mag(T.internalField() - T.prevIter().internalField()));
    const scalar residualT =
        absResidualT
       /max
        (
            gMax(mag(T.internalField() - T.oldTime().internalField())),
            SMALL
        );

    const scalar residualD =
        gMax
        (
            mag(D.internalField() - D.prevIter().internalField())
           /max
            (
                gMax(mag(D.internalField() - D.oldTime().internalField())),
                SMALL
            )
        );

    // Calculate material residual
    const scalar materialResidual = mechanical().residual();

    // If one of the residuals has converged to an order of magnitude
    // less than the tolerance then consider the solution converged
    // force at leaast 1 outer iteration and the material law must be converged
    if (iCorr > 1 && materialResidual < materialTol())
    {
        bool convergedD = false;
        bool convergedT = false;

        if
        (
            (
                solverPerfD.initialResidual() < solutionTol()
             && residualD < solutionTol()
            )
         || solverPerfD.initialResidual() < alternativeTol()
         || residualD < alternativeTol()
        )
        {
            convergedD = true;
        }

        if
        (
            (
                solverPerfT.initialResidual() < solutionTol()
             && residualT < solutionTol()
            )
         || solverPerfT.initialResidual() < alternativeTol()
         || residualT < alternativeTol()
         || absResidualT < absTTol_
        )
        {
            convergedT = true;
        }


        if (convergedD && convergedT)
        {
            Info<< "    The residuals have converged" << endl;
            converged = true;
        }
    }

    // Print residual information
    if (iCorr == 0)
    {
        Info<< "    Corr, res (T & D), relRes (T & D), matRes, iters (T & D), enforceLinear"
            << endl;
    }
    else if (iCorr % infoFrequency() == 0 || converged)
    {
        Info<< "    " << iCorr
            << ", " << solverPerfT.initialResidual()
            << ", " << solverPerfD.initialResidual()
            << ", " << residualT
            << ", " << residualD
            << ", " << materialResidual
            << ", " << solverPerfT.nIterations()
            << ", " << solverPerfD.nIterations()
            << ", " << enforceLinear() << endl;

        if (converged)
        {
            Info<< endl;
        }
    }
    else if (iCorr == nCorr() - 1)
    {
        maxIterReached()++;
        Warning
            << "Max iterations reached within the enery-momentum loop" << endl;
    }

    return converged;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

unsThermalNonLinGeomTotalLagSolid::unsThermalNonLinGeomTotalLagSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    sigmaf_
    (
        IOobject
        (
            "sigmaf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    ),
    gradDf_
    (
        IOobject
        (
            "grad(" + D().name() + ")f",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimless, tensor::zero)
    ),
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
    Ff_
    (
        IOobject
        (
            "Ff",
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
    Finvf_
    (
        IOobject
        (
            "Finvf",
            runTime.timeName(),
            mesh(),
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
            runTime.timeName(),
            mesh(),
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
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        det(Ff_)
    ),
    thermal_(mesh()),
    rhoC_
    (
        IOobject
        (
            "rhoC",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        thermal_.C()*mechanical().rho()
    ),
    kappa_(thermal_.k()),
    T_
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
    gradT_
    (
        IOobject
        (
            "grad(T)",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimTemperature/dimLength, vector::zero)
    ),
    absTTol_
    (
        solidModelDict().lookupOrDefault<scalar>
        (
            "absoluteTemperatureTolerance",
            1e-06
        )
    ),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_),
    nonLinear_(solidModelDict().lookupOrDefault<Switch>("nonLinear", true)),
    debug_(solidModelDict().lookupOrDefault<Switch>("debug", false)),
    K_
    (
        solidModelDict().lookupOrDefault<dimensionedScalar>
        (
            "K",
            dimensionedScalar("K", dimless/dimTime, 0)
        )
    ),
    relativeTol_
    (
        solidModelDict().lookupOrDefault<scalar>
        (
            "solutionTolerance",
            solutionTol()
        )
    )
{
    DisRequired();

    // Store T old time
    T_.oldTime();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool unsThermalNonLinGeomTotalLagSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    int iCorr = 0;
    lduSolverPerformance solverPerfD;
    lduSolverPerformance solverPerfT;
    blockLduMatrix::debug = 0;

    // Reset enforceLinear switch
    enforceLinear() = false;

    Info<< "Solving coupled energy and displacements equation for T and D"
        << endl;

    // Momentum-energy coupling outer loop
    do
    {
        if (blockLduMatrix::debug)
        {
            Info<< "Time: " << runTime().timeName()
                << ", outer iteration: " << iCorr << endl;
        }

        // Store fields for under-relaxation and residual calculation
        T().storePrevIter();

        // Heat equation
        fvScalarMatrix TEqn
        (
            rhoC_*fvm::ddt(T_)
         == fvm::laplacian(kappa_, T_, "laplacian(k,T)")
          + (sigma() && fvc::grad(U()))
        );

        // Under-relaxation the linear system
        TEqn.relax();

        // Solve the linear system
        solverPerfT = TEqn.solve();

        // Under-relax the field
        T_.relax();

        // Update gradient of temperature
        gradT_ = fvc::grad(T_);

        // Store previous iteration to allow under-relaxation and residual
        // calculation
        D().storePrevIter();

        // Construct momentum equation in total Lagrangian form where gradients
        // are calculated directly at the faces
        fvVectorMatrix DEqn
        (
            rho()*fvm::d2dt2(D())
         == fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
          - fvc::laplacian(impKf_, D(), "laplacian(DD,D)")
          + fvc::div((Jf_*Finvf_.T() & mesh().Sf()) & sigmaf_)
          + rho()*g()
        );

        // Add damping
        if (K_.value() > SMALL)
        {
            DEqn += K_*rho()*fvm::ddt(D());
        }

        // Enforce linear to improve convergence
        if (enforceLinear())
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
        solverPerfD = DEqn.solve();

        // Under-relax displacement field
        relaxField(D(), iCorr);

        // Interpolate D to pointD
        mechanical().interpolate(D(), pointD(), false);

        // Update gradient of displacement
        mechanical().grad(D(), pointD(), gradD());
        mechanical().grad(D(), pointD(), gradDf_);

        // Update gradient of displacement increment
        gradDD() = gradD() - gradD().oldTime();

        // Total deformation gradient
        Ff_ = I + gradDf_.T();

        // Inverse of the deformation gradient
        Finvf_ = inv(Ff_);

        // Jacobian of the deformation gradient
        Jf_ = det(Ff_);

        // Check if outer loops are diverging
        if (nonLinear_ && !enforceLinear())
        {
            checkEnforceLinear(Jf_);
        }

        // Calculate the stress using run-time selectable mechanical law
        mechanical().correct(sigmaf_);
    }
    while
    (
        !converged(iCorr, solverPerfD, solverPerfT, D(), T_)
     && ++iCorr < nCorr()
    );

    // Velocity
    U() = fvc::ddt(D());

    // Total deformation gradient
    F_ = I + gradD().T();

    // Inverse of the deformation gradient
    Finv_ = inv(F_);

    // Jacobian of the deformation gradient
    J_ = det(F_);

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());

    // Increment of displacement
    DD() = D() - D().oldTime();

    // Increment of point displacement
    pointDD() = pointD() - pointD().oldTime();

    blockLduMatrix::debug = 1;

    if (nonLinear_ && enforceLinear())
    {
        return false;
    }

    return true;
}


tmp<vectorField> unsThermalNonLinGeomTotalLagSolid::tractionBoundarySnGrad
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

    // Patch unit normals (initial configuration)
    const vectorField n = patch.nf();

    if (enforceLinear())
    {
        // Return patch snGrad
        return tmp<vectorField>
        (
            new vectorField
            (
                (
                    (traction - n*pressure)
                  - (n & sigma)
                  + (n & (impK*gradD))
                )*rImpK
            )
        );
    }
    else
    {
        // Patch total deformation gradient inverse
        const tensorField& Finv = Finvf_.boundaryField()[patchID];

        // Patch total Jacobian
        const scalarField& J = Jf_.boundaryField()[patchID];

        // Patch unit normals (deformed configuration)
        const vectorField nCurrent = J*Finv.T() & n;

        // Return patch snGrad
        return tmp<vectorField>
        (
            new vectorField
            (
                (
                    (traction - nCurrent*pressure)
                  - (nCurrent & sigma)
                  + (n & (impK*gradD))
                )*rImpK
            )
        );
    }
}


void unsThermalNonLinGeomTotalLagSolid::writeFields(const Time& runTime)
{
    Info<< "Max T = " << max(T_).value() << nl
        << "Min T = " << min(T_).value() << endl;

    // Heat flux
    volVectorField heatFlux
    (
        IOobject
        (
            "heatFlux",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
       -kappa_*gradT_
    );

    Info<< "Max magnitude of heat flux = " << max(mag(heatFlux)).value()
        << endl;

    solidModel::writeFields(runTime);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
