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

#include "thermalLinGeomSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "mixedFvPatchFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(thermalLinGeomSolid, 0);
addToRunTimeSelectionTable(physicsModel, thermalLinGeomSolid, solid);
addToRunTimeSelectionTable(solidModel, thermalLinGeomSolid, dictionary);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

bool thermalLinGeomSolid::converged
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
        Info<< "    Corr, res (T & D), relRes (T & D), matRes, iters (T & D)"
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
            << ", " << solverPerfD.nIterations() << endl;

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

thermalLinGeomSolid::thermalLinGeomSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
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
    DiffusionNo_(0)
{
    DisRequired();
    TisRequired();

    // Store T old time
    T().oldTime();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar& thermalLinGeomSolid::DiffusionNo()
{
    //- calculate solid Diffusion number
    DiffusionNo_ = 0.0;
    scalar meanDiffusionNo = 0.0;

    //- Can have fluid domains with 0 cells so do not test.
    if (mesh().nInternalFaces())
    {
           surfaceScalarField kRhoCbyDelta =
               mesh().surfaceInterpolation::deltaCoeffs()
             * fvc::interpolate(kappa_)
             / fvc::interpolate(rhoC_);

           DiffusionNo_ = max(kRhoCbyDelta.internalField())*runTime().deltaT().value();

           meanDiffusionNo = (average(kRhoCbyDelta)).value()*runTime().deltaT().value();
    }

    Info<< "Diffusion Number mean: " << meanDiffusionNo
        << " max: " << DiffusionNo_ << endl;

    return DiffusionNo_;
}


tmp<scalarField> thermalLinGeomSolid::patchHeatFlux
(
    const label patchID
) const
{
    tmp<scalarField> ttF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    ttF() = kappa_.boundaryField()[patchID]
          * T().boundaryField()[patchID].snGrad();

    return ttF;
}


tmp<scalarField> thermalLinGeomSolid::patchTemperature
(
    const label patchID
) const
{
    tmp<scalarField> tT
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    tT() = T().boundaryField()[patchID].patchInternalField();

    return tT;
}


tmp<scalarField> thermalLinGeomSolid::patchKappaDelta
(
    const label patchID
) const
{
    tmp<scalarField> tKD
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    tKD() = kappa_.boundaryField()[patchID]
          * mesh().boundary()[patchID].deltaCoeffs();

    return tKD;
}


void thermalLinGeomSolid::setTemperature
(
    const label patchID,
    const scalarField& nbrFaceZoneTemperature,
    const scalarField& nbrFaceZoneKappaDelta
)
{
    if
    (
        T().boundaryField()[patchID].type()
     != mixedFvPatchScalarField::typeName
    )
    {
        FatalErrorIn
        (
            "void thermalLinGeomSolid::setTemperature\n"
            "(\n"
            "    const label,\n"
            "    const scalarField&,\n"
            "    const scalarField&\n"
            ")"
        )
            << "Bounary condition on " << T().name()
                <<  " is "
                << T().boundaryField()[patchID].type()
                << " for patch " << mesh().boundary()[patchID].name()
                << ", instead of "
                << mixedFvPatchScalarField::typeName
                << abort(FatalError);
    }

    scalarField nbrPatchTemperature =
	globalPatch().globalFaceToPatch(nbrFaceZoneTemperature);

    scalarField nbrPatchKappaDelta =
	globalPatch().globalFaceToPatch(nbrFaceZoneKappaDelta);

    mixedFvPatchScalarField& patchT =
        refCast<mixedFvPatchScalarField>
        (
            T().boundaryField()[patchID]
        );

    patchT.refValue() = nbrPatchTemperature;

    patchT.refGrad() = 0.0;

    patchT.valueFraction() =
        nbrPatchKappaDelta
      / (nbrPatchKappaDelta + patchKappaDelta(patchID));

    patchT.updateCoeffs();
}


bool thermalLinGeomSolid::evolve()
{
    Info<< "Evolving thermal solid solver" << endl;

    int iCorr = 0;
    lduSolverPerformance solverPerfD;
    lduSolverPerformance solverPerfT;
    blockLduMatrix::debug = 0;

    // calculate Diffusion number
    DiffusionNo();

    Info<< "Solving coupled energy and displacements equation for T and D"
        << endl;

    // Momentum-energy coupling outer loop
    do
    {
        // Store fields for under-relaxation and residual calculation
        T().storePrevIter();

        // Heat equation
        fvScalarMatrix TEqn
        (
            rhoC_*fvm::ddt(T())
         == fvm::laplacian(kappa_, T(), "laplacian(k,T)")
          + (sigma() && fvc::grad(U()))
          - fvm::SuSp(-thermal_.S()/T(), T())
        );

        // Under-relaxation the linear system
        TEqn.relax();

        // Solve the linear system
        solverPerfT = TEqn.solve();

        // Under-relax the field
        T().relax();

        // Update gradient of temperature
        gradT() = fvc::grad(T());

        // Store fields for under-relaxation and residual calculation
        D().storePrevIter();

        // Linear momentum equation total displacement form
        fvVectorMatrix DEqn
        (
            rho()*fvm::d2dt2(D())
         == fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
          - fvc::laplacian(impKf_, D(), "laplacian(DD,D)")
          + fvc::div(sigma(), "div(sigma)")
          + rho()*g()
          + mechanical().RhieChowCorrection(D(), gradD())
        );

        // Under-relaxation the linear system
        DEqn.relax();

        // Enforce any cell displacements
        solidModel::setCellDisps(DEqn);

        // Solve the linear system
        solverPerfD = DEqn.solve();

        // Under-relax the field
        relaxField(D(), iCorr);

        // Update increment of displacement
        DD() = D() - D().oldTime();

        // Update velocity
        U() = fvc::ddt(D());

        // Update gradient of displacement
        mechanical().grad(D(), gradD());

        // Update gradient of displacement increment
        gradDD() = gradD() - gradD().oldTime();

        // Calculate the stress using run-time selectable mechanical law
        mechanical().correct(sigma());

        // Update impKf to improve convergence
        // Note: impK and rImpK are not updated as they are used for traction
        // boundaries
        if (iCorr % 10 == 0)
        {
            impKf_ = mechanical().impKf();
        }
    }
    while
    (
        !converged(iCorr, solverPerfD, solverPerfT, D(), T())
     && ++iCorr < nCorr()
    );

    // Interpolate cell displacements to vertices
    mechanical().interpolate(D(), pointD());

    // Increment of displacement
    DD() = D() - D().oldTime();

    // Increment of point displacement
    pointDD() = pointD() - pointD().oldTime();

    blockLduMatrix::debug = 1;

    return true;
}


tmp<vectorField> thermalLinGeomSolid::tractionBoundarySnGrad
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


void thermalLinGeomSolid::writeFields(const Time& runTime)
{
    Info<< "Max T = " << max(T()).value() << nl
        << "Min T = " << min(T()).value() << endl;

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
       -kappa_*gradT()
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
