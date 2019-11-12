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

#include "unsWeakThermalNonLinGeomTotalLagSolid.H"
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

defineTypeNameAndDebug(unsWeakThermalNonLinGeomTotalLagSolid, 0);
addToRunTimeSelectionTable
(
    physicsModel, unsWeakThermalNonLinGeomTotalLagSolid, solid
);
addToRunTimeSelectionTable
(
    solidModel, unsWeakThermalNonLinGeomTotalLagSolid, dictionary
);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

// Note: Sep 18 2019
// The total Lagrangian implementation does not move the mesh
void unsWeakThermalNonLinGeomTotalLagSolid::DiffusionNo()
{
    //- calculate solid Diffusion number
    DiffusionNum() = 0.0;
    scalar meanDiffusionNum = 0.0;

    //- Can have fluid domains with 0 cells so do not test.
    if (mesh().nInternalFaces())
    {
           surfaceScalarField kRhoCbyDelta =
               mesh().surfaceInterpolation::deltaCoeffs()
             * fvc::interpolate(kappa_)
             / fvc::interpolate(rhoC_);

           DiffusionNum() =
               max(kRhoCbyDelta.internalField())*runTime().deltaT().value();

           meanDiffusionNum =
               (average(kRhoCbyDelta)).value()*runTime().deltaT().value();
    }

    Info<< "Diffusion Number mean: " << meanDiffusionNum
        << " max: " << DiffusionNum() << endl;
}


bool unsWeakThermalNonLinGeomTotalLagSolid::converged
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


const standAlonePatch& unsWeakThermalNonLinGeomTotalLagSolid::boundaryPatchCurrent
(
    const label patchID
) const
{
    deformedPatchPtr_ =
        new standAlonePatch
        (
            mesh().boundaryMesh()[patchID].localFaces(),
            mesh().boundaryMesh()[patchID].localPoints()
        );
    standAlonePatch& deformedPatch = *deformedPatchPtr_;

    // Calculate the deformed points
    const pointField deformedPoints =
        mechanical().volToPoint().interpolate
        (
            mesh().boundaryMesh()[patchID],
            D()
        )
      + mesh().boundaryMesh()[patchID].localPoints();

    // Move the standAlonePatch points
    const_cast<pointField&>(deformedPatch.points()) = deformedPoints;

    // Return the boundary patch in its deformed configuration
    return *deformedPatchPtr_;
}


tmp<scalarField> unsWeakThermalNonLinGeomTotalLagSolid::deltaCoeffsCurrent
(
    const label patchID
) const
{
    tmp<scalarField> tdeltaCoeffsCurrent
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    // Patch delta vector (initial configuration)
    const vectorField delta = mesh().boundary()[patchID].delta();

    // Patch delta vector (deformed configuration)
    const vectorField deltaCurrent =
    (
        D().boundaryField()[patchID]
      - D().boundaryField()[patchID].patchInternalField()
    ) + delta;

    // Note: we have two options to calculate deformed unit normal
    // 1. Use the cell-centred deformation gradient field 'F'.
    // 2. Calculate the deformed normals by interpolating displacements
    //    to the points and calculating the normals on the deformed patch.
    // For now we use method 1. as the result of the method 2. is not consistent,
    // perhaps due the governing equation formulation

    // Patch unit normals (initial configuration)
    const vectorField n = mesh().boundary()[patchID].nf();

    // Patch total deformation gradient inverse
    const tensorField& FinvBf = Finvf().boundaryField()[patchID];

    // Patch total Jacobian
    const scalarField& JBf = Jf().boundaryField()[patchID];

    // Patch unit normals (deformed configuration)
    const vectorField nCurrent = JBf*FinvBf.T() & n;

    // Patch unit normals (deformed configuration)
    // const vectorField& nCurrent
    // (
    //     boundaryPatchCurrent(patchID).faceNormals()
    // );

    forAll(tdeltaCoeffsCurrent(), faceI)
    {
        tdeltaCoeffsCurrent()[faceI] =
            scalar(1) / max
            (
                nCurrent[faceI] & deltaCurrent[faceI],
                0.05*mag(deltaCurrent[faceI])
            );
    }

    return tdeltaCoeffsCurrent;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

unsWeakThermalNonLinGeomTotalLagSolid::unsWeakThermalNonLinGeomTotalLagSolid
(
    Time& runTime,
    const word& region
)
:
    unsNonLinGeomTotalLagSolid(runTime, region),
    deformedPatchPtr_(NULL),
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
    )
{
    TisRequired();

    // Store T old time
    T().oldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

unsWeakThermalNonLinGeomTotalLagSolid::~unsWeakThermalNonLinGeomTotalLagSolid()
{
    deleteDemandDrivenData(deformedPatchPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> unsWeakThermalNonLinGeomTotalLagSolid::patchHeatFlux
(
    const label patchID
) const
{
    tmp<scalarField> thF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    thF() =
        kappa_.boundaryField()[patchID]
      * T().boundaryField()[patchID].snGrad();

    return thF;
}


tmp<scalarField> unsWeakThermalNonLinGeomTotalLagSolid::patchTemperature
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


tmp<scalarField> unsWeakThermalNonLinGeomTotalLagSolid::patchKappaDelta
(
    const label patchID
) const
{
    tmp<scalarField> tKD
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    tKD() =
        kappa_.boundaryField()[patchID]*deltaCoeffsCurrent(patchID);

    return tKD;
}


bool unsWeakThermalNonLinGeomTotalLagSolid::evolve()
{
    Info<< "Evolving thermal solid solver" << endl;

    int iCorr = 0;
    lduSolverPerformance solverPerfT;
    blockLduMatrix::debug = 0;

    Info<< "Solving the energy equation for T" << endl;

    // Energy equation non-orthogonality corrector loop
    do
    {
        // Store fields for under-relaxation and residual calculation
        T().storePrevIter();

        // Heat equation
        fvScalarMatrix TEqn
        (
            rhoC_*fvm::ddt(T())
         == fvm::laplacian(kappa_, T(), "laplacian(k,T)")
          - fvc::laplacian(kappa_, T(), "laplacian(k,T)")
          + fvc::div
            (
                (Jf()*Finvf().T() & mesh().Sf())
              & fvc::interpolate(kappa_*gradT())
            )
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

        // Update rho*C
        rhoC_ = rho()*thermal_.C();
    }
    while (!converged(iCorr, solverPerfT, T()) && ++iCorr < nCorr());

    // Now solve the momentum equation
    // Note: a thermal elastic mechanical law should be chosen if the effect of
    // thermal stress is to be included in the momentum equation.
    unsNonLinGeomTotalLagSolid::evolve();

    Info<< "Solid temperature min/max(T) = " << min(T()).value()
	<< ", " << max(T()).value() << " [K]" << endl;

    blockLduMatrix::debug = 1;

    return true;
}


void unsWeakThermalNonLinGeomTotalLagSolid::writeFields(const Time& runTime)
{
    Info<< "Max T = " << max(T()).value() << nl
        << "Min T = " << min(T()).value() << endl;

    volScalarField wallHeatFlux
    (
        IOobject
        (
            "wallHeatFlux",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimMass/pow3(dimTime), 0.0)
    );

    volScalarField::GeometricBoundaryField& wallHeatFluxBf =
        wallHeatFlux.boundaryField();

    Info<< "\nWall heat fluxes [W]" << endl;
    forAll(wallHeatFluxBf, patchI)
    {
        if (mesh().boundary()[patchI].isWall())
        {
            wallHeatFluxBf[patchI] =
                kappa_.boundaryField()[patchI]
              * T().boundaryField()[patchI].snGrad();

            Info<< mesh().boundary()[patchI].name()
                << " = "
                << gSum
                   (
                       wallHeatFluxBf[patchI]
                     * mag(patchCurrentSf(patchI))
                   )
                << endl;
        }
    }
    Info<< endl;

    unsNonLinGeomTotalLagSolid::writeFields(runTime);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
