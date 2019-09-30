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

#include "thermalSolid.H"
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

defineTypeNameAndDebug(thermalSolid, 0);
addToRunTimeSelectionTable(physicsModel, thermalSolid, solid);
addToRunTimeSelectionTable(solidModel, thermalSolid, dictionary);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

bool thermalSolid::converged
(
    const int iCorr,
    const lduSolverPerformance& solverPerfT,
    const volScalarField& T
)
{
    // We will check a number of different residuals for convergence
    bool converged = false;

    // Calculate relative residuals
    const scalar absResidualT =
        gMax(mag(T.internalField() - T.oldTime().internalField()));
    const scalar residualT =
        gMax
        (
            mag(T.internalField() - T.prevIter().internalField())
           /max
            (
                gMax(mag(T.internalField() - T.oldTime().internalField())),
                SMALL
            )
        );

    // If one of the residuals has converged to an order of magnitude
    // less than the tolerance then consider the solution converged
    // force at leaast 1 outer iteration and the material law must be converged
    if (iCorr > 1)
    {
        if (absResidualT < absTTol_)
        {
            Info<< "    T has converged to within the " << absTTol_
                << " degrees" << endl;
            converged = true;
        }
        else if
        (
            solverPerfT.initialResidual() < solutionTol()
         && residualT < solutionTol()
        )
        {
            Info<< "    Both residuals have converged" << endl;
            converged = true;
        }
        else if
        (
            residualT < alternativeTol()
        )
        {
            Info<< "    The relative residual has converged" << endl;
            converged = true;
        }
        else if
        (
            solverPerfT.initialResidual() < alternativeTol()
        )
        {
            Info<< "    The solver residual has converged" << endl;
            converged = true;
        }
        else
        {
            converged = false;
        }
    }

    // Print residual information
    if (iCorr == 0)
    {
        Info<< "    Corr, res, relRes, iters" << endl;
    }
    else if (iCorr % infoFrequency() == 0 || converged)
    {
        Info<< "    " << iCorr
            << ", " << solverPerfT.initialResidual()
            << ", " << residualT
            << ", " << solverPerfT.nIterations() << endl;

        if (converged)
        {
            Info<< endl;
        }
    }
    else if (iCorr == nCorr() - 1)
    {
        maxIterReached()++;
        Warning
            << "Max iterations reached within temperature loop" << endl;
    }

    return converged;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermalSolid::thermalSolid
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
    DiffusionNo_(0)
{
    TisRequired();

    // Store T old time
    T().oldTime();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar& thermalSolid::DiffusionNo()
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


tmp<scalarField> thermalSolid::patchHeatFlux
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


tmp<scalarField> thermalSolid::patchTemperature
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


tmp<scalarField> thermalSolid::patchKappaDelta
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


void thermalSolid::setTemperature
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
            "void thermalSolid::setTemperature\n"
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


bool thermalSolid::evolve()
{
    Info<< "Evolving thermal solid solver" << endl;

    int iCorr = 0;
    lduSolverPerformance solverPerfT;
    blockLduMatrix::debug = 0;

    // calculate Diffusion number
    DiffusionNo();

    // energy equation outer loop
    do
    {
        // Store fields for under-relaxation and residual calculation
        T().storePrevIter();

        // Heat equation
        fvScalarMatrix TEqn
        (
            rhoC_*fvm::ddt(T())
          - fvm::laplacian(kappa_, T(), "laplacian(k,T)")
          + fvm::SuSp(-thermal_.S()/T(), T())
        );

        // Under-relaxation the linear system
        TEqn.relax();

        // Solve the linear system
        solverPerfT = TEqn.solve();

        // Under-relax the field
        T().relax();

        // Update gradient of temperature
        gradT() = fvc::grad(T());
    }
    while
    (
        !converged(iCorr, solverPerfT, T())
     && ++iCorr < nCorr()
    );

    blockLduMatrix::debug = 1;

    return true;
}


tmp<vectorField> thermalSolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    notImplemented
    (
        "thermalSolid::tractionBoundarySnGrad\n"
        "(\n"
        "    const vectorField&,\n"
        "    const scalarField&,\n"
        "    const fvPatch&\n"
        ") const"
    );

    return tmp<vectorField>
    (
        new vectorField(patch.size(), vector::zero)
    );
}


void thermalSolid::writeFields(const Time& runTime)
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

    physicsModel::writeFields(runTime);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
