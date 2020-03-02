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

void thermalSolid::solidRegionDiffNo()
{
    surfaceScalarField kapparhoCpbyDelta
    (
        sqr(mesh().surfaceInterpolation::deltaCoeffs())
       *fvc::interpolate(kappa_)
       /fvc::interpolate(C_*rho_)
    );

    const scalar DiNum =
        max(kapparhoCpbyDelta).value()*runTime().deltaTValue();

    const scalar meanDiNum =
        average(kapparhoCpbyDelta).value()*runTime().deltaTValue();

    Info<< "Diffusion Number mean: " << meanDiNum
        << " max: " << DiNum << endl;
}


bool thermalSolid::converged
(
    const int iCorr,
#ifdef OPENFOAMESIORFOUNDATION
    const SolverPerformance<scalar>& solverPerfT,
#else
    const lduSolverPerformance& solverPerfT,
#endif
    const volScalarField& T
)
{
    // We will check a number of different residuals for convergence
    bool converged = false;

    // Calculate relative residuals
    const scalar absResidualT =
#ifdef OPENFOAMESIORFOUNDATION
        gMax(mag(T.primitiveField() - T.oldTime().primitiveField()));
#else
        gMax(mag(T.internalField() - T.oldTime().internalField()));
#endif
    const scalar residualT =
        gMax
        (
#ifdef OPENFOAMESIORFOUNDATION
            mag(T.primitiveField() - T.prevIter().primitiveField())
           /max
            (
                gMax(mag(T.primitiveField() - T.oldTime().primitiveField())),
                SMALL
            )
#else
            mag(T.internalField() - T.prevIter().internalField())
           /max
            (
                gMax(mag(T.internalField() - T.oldTime().internalField())),
                SMALL
            )
#endif
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
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimTemperature/dimLength, vector::zero)
    ),
    rho_(mechanical().rho()),
    C_(thermal_.C()),
    kappa_(thermal_.k()),
    absTTol_
    (
        solidModelDict().lookupOrDefault<scalar>
        (
            "absoluteTemperatureTolerance",
            1e-06
        )
    )
{
    // Store T old time
    T_.oldTime();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> thermalSolid::patchHeatFlux
(
    const label patchID
) const
{
    tmp<scalarField> thF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

#ifdef OPENFOAMESIORFOUNDATION
    thF.ref() =
#else
    thF() =
#endif
        kappa_.boundaryField()[patchID]
      * T_.boundaryField()[patchID].snGrad();

    return thF;
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

#ifdef OPENFOAMESIORFOUNDATION
    tT.ref() =
#else
    tT() =
#endif
        T_.boundaryField()[patchID].patchInternalField();

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

#ifdef OPENFOAMESIORFOUNDATION
    tKD.ref() =
#else
    tKD() =
#endif
        kappa_.boundaryField()[patchID]
      * mesh().boundary()[patchID].deltaCoeffs();

    return tKD;
}


bool thermalSolid::evolve()
{
    Info<< "Evolving thermal solid solver" << endl;

    int iCorr = 0;
#ifdef OPENFOAMESIORFOUNDATION
    SolverPerformance<scalar> solverPerfT;
    SolverPerformance<scalar>::debug = 0;
#else
    lduSolverPerformance solverPerfT;
    blockLduMatrix::debug = 0;
#endif

    solidRegionDiffNo();

    // energy equation outer loop
    do
    {
        // Store fields for under-relaxation and residual calculation
        T_.storePrevIter();

        // Heat equation
        fvScalarMatrix TEqn
        (
            fvm::ddt(rho_*C_, T_)
          - fvm::laplacian(kappa_, T_, "laplacian(kappa,T)")
        );

        // Under-relaxation the linear system
        TEqn.relax();

        // Solve the linear system
        solverPerfT = TEqn.solve();

        // Under-relax the field
        T_.relax();

        // Update gradient of temperature
        gradT_ = fvc::grad(T_);
    }
    while
    (
        !converged(iCorr, solverPerfT, T_)
     && ++iCorr < nCorr()
    );

    Info<< "Solid temperature min/max(T) = " << min(T_).value()
	<< ", " << max(T_).value() << " [K]" << endl;

#ifndef OPENFOAMESIORFOUNDATION
    blockLduMatrix::debug = 1;
#endif

    return true;
}


void thermalSolid::writeFields(const Time& runTime)
{
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

    Info<< "Max magnitude of heat flux = "
        << max(mag(heatFlux)).value() << endl;

    physicsModel::writeFields(runTime);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
