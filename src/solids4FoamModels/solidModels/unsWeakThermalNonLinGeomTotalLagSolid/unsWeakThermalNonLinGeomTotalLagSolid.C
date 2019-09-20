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
#include "mixedFvPatchFields.H"


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


tmp<vectorField> unsWeakThermalNonLinGeomTotalLagSolid::currentFaceNormal
(
    const label patchID
) const
{
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

    standAlonePatch deformedPatch =
        standAlonePatch
        (
            mesh().boundaryMesh()[patchID].localFaces(),
            mesh().boundaryMesh()[patchID].localPoints()
        );

    // Calculate the deformed points
    const pointField deformedPoints =
        mechanical().volToPoint().interpolate
        (
            mesh().boundaryMesh()[patchID],
            DD()
        )
      + mesh().boundaryMesh()[patchID].localPoints();

    // Move the standAlonePatch points
    const_cast<pointField&>(deformedPatch.points()) = deformedPoints;

    // Patch unit normals (deformed configuration)
    return deformedPatch.faceNormals();
}


tmp<scalarField> unsWeakThermalNonLinGeomTotalLagSolid::currentDeltaCoeffs
(
    const label patchID
) const
{
    tmp<scalarField> tcurrentDelta
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    // Patch delta vector (initial configuration)
    const vectorField delta = mesh().boundary()[patchID].delta();

    // Patch delta vector (deformed configuration)
    const vectorField deltaCurrent =
    (
        DD().boundaryField()[patchID]
      - DD().boundaryField()[patchID].patchInternalField()
    ) + delta;

    // Patch unit normals (deformed configuration)
    const vectorField& nCurrent = currentFaceNormal(patchID);

    forAll(tcurrentDelta(), faceI)
    {
        tcurrentDelta()[faceI] = max
            (
                nCurrent[faceI] & deltaCurrent[faceI],
                0.05*mag(deltaCurrent[faceI])
            );
    }

    return scalar(1.0)/tcurrentDelta;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

unsWeakThermalNonLinGeomTotalLagSolid::unsWeakThermalNonLinGeomTotalLagSolid
(
    Time& runTime,
    const word& region
)
:
    unsNonLinGeomTotalLagSolid(runTime, region),
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

// Needs correction: Sep 18 2019
// The total Lagrangian implementation does not move the mesh
scalar& unsWeakThermalNonLinGeomTotalLagSolid::DiffusionNo()
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


tmp<scalarField> unsWeakThermalNonLinGeomTotalLagSolid::patchThermalFlux
(
    const label patchID
) const
{
    tmp<scalarField> ttF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    // corrected snGrad (deformed configuration)
    const scalarField snGradT =
    (
        T().boundaryField()[patchID]
      - T().boundaryField()[patchID].patchInternalField()
    ) * currentDeltaCoeffs(patchID);

    ttF() = kappa_.boundaryField()[patchID]*snGradT;

    return ttF;
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


tmp<scalarField> unsWeakThermalNonLinGeomTotalLagSolid::patchKDelta
(
    const label patchID
) const
{
    tmp<scalarField> tKD
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    tKD() = kappa_.boundaryField()[patchID]*currentDeltaCoeffs(patchID);

    return tKD;
}


void unsWeakThermalNonLinGeomTotalLagSolid::setTemperature
(
    const label patchID,
    const scalarField& nbrFaceZoneTemperature,
    const scalarField& nbrFaceZoneKDelta
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
            "void unsWeakThermalNonLinGeomTotalLagSolid::setTemperature\n"
            "(\n"
            "    const label,\n"
            "    const scalarField&,\n"
            "    const scalarField&\n"
            ")"
        )
            << "Bounary condition on " << T().name()
                <<  " is "
                << T().boundaryField()[patchID].type()
                << "for patch" << mesh().boundary()[patchID].name()
                << ", instead of "
                << mixedFvPatchScalarField::typeName
                << abort(FatalError);
    }

    scalarField nbrPatchTemperature =
	globalPatch().globalFaceToPatch(nbrFaceZoneTemperature);

    scalarField nbrPatchKDelta =
	globalPatch().globalFaceToPatch(nbrFaceZoneKDelta);

    mixedFvPatchScalarField& patchT =
        refCast<mixedFvPatchScalarField>
        (
            T().boundaryField()[patchID]
        );

    patchT.refValue() = nbrPatchTemperature;
    patchT.refGrad() = 0.0;
    patchT.valueFraction() = nbrPatchKDelta / (nbrPatchKDelta + patchKDelta(patchID));
    patchT.updateCoeffs();
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
          + fvc::div(J()*kappa_*gradT() & Finv().T(), "div(k*grad(T))")
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
    }
    while (!converged(iCorr, solverPerfT, T()) && ++iCorr < nCorr());

    // Now solve the momentum equation
    // Note: a thermal elastic mechanical law should be chosen if the effect of
    // thermal stress is to be included in the momentum equation.
    unsNonLinGeomTotalLagSolid::evolve();

    blockLduMatrix::debug = 1;

    return true;
}


void unsWeakThermalNonLinGeomTotalLagSolid::writeFields(const Time& runTime)
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

    unsNonLinGeomTotalLagSolid::writeFields(runTime);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
