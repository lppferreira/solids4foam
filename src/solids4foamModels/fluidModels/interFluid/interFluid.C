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

#include "interFluid.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "findRefCell.H"
#include "adjustPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(interFluid, 0);
addToRunTimeSelectionTable(physicsModel, interFluid, fluid);
addToRunTimeSelectionTable(fluidModel, interFluid, dictionary);


// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

void interFluid::solveAlphaEqnSubCycle()
{
    dictionary pimple = mesh().solutionDict().subDict("PIMPLE");

    const label nAlphaCorr
    (
        readLabel(pimple.lookup("nAlphaCorr"))
    );

    const label nAlphaSubCycles
    (
        readLabel(pimple.lookup("nAlphaSubCycles"))
    );

    if (nAlphaSubCycles > 1)
    {
        const dimensionedScalar totalDeltaT = runTime().deltaT();
        surfaceScalarField rhoPhiSum = 0.0*rhoPhi_;

        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha1_, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            solveAlphaEqn(nAlphaCorr);
            rhoPhiSum += (runTime().deltaT()/totalDeltaT)*rhoPhi_;
        }

        rhoPhi_ = rhoPhiSum;
    }
    else
    {
        solveAlphaEqn(nAlphaCorr);
    }

    interface_.correct();

    rho_ == alpha1_*rho1_ + (scalar(1) - alpha1_)*rho2_;
}


void interFluid::solveAlphaEqn(const label nAlphaCorr)
{
    const word alphaScheme("div(phi,alpha)");
    const word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phic = mag(phi()/mesh().magSf());
    phic = min(interface_.cAlpha()*phic, max(phic));
    surfaceScalarField phir = phic*interface_.nHatf();

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        surfaceScalarField phiAlpha =
            fvc::flux
            (
                phi(),
                alpha1_,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, scalar(1) - alpha1_, alpharScheme),
                alpha1_,
                alpharScheme
            );

        MULES::explicitSolve(alpha1_, phi(), phiAlpha, 1, 0);

        rhoPhi_ = phiAlpha*(rho1_ - rho2_) + phi()*rho2_;
    }

    Info<< "Liquid phase volume fraction = "
        << alpha1_.weightedAverage(mesh().V()).value()
        << "  Min(alpha1) = " << min(alpha1_).value()
        << "  Max(alpha1) = " << max(alpha1_).value()
        << endl;
}


tmp<fvVectorMatrix> interFluid::solveUEqn()
{
    dictionary pimple = mesh().solutionDict().subDict("PIMPLE");

    const surfaceScalarField muEff
    (
        "muEff",
        twoPhaseProperties_.muf()
      + fvc::interpolate(rho_*turbulence_->nut())
    );

    tmp<fvVectorMatrix> tUEqn
    (
        new fvVectorMatrix
        (
            fvm::ddt(rho_, U())
          + fvm::div(rhoPhi_, U())
          - fvm::laplacian(muEff, U())
          - (fvc::grad(U()) & fvc::grad(muEff))
        )
    );
    fvVectorMatrix& UEqn = tUEqn();

    UEqn.relax();

    if (pimple.lookupOrDefault<Switch>("momentumPredictor", true))
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    fvc::interpolate(interface_.sigmaK())*fvc::snGrad(alpha1_)
                  - ghf_*fvc::snGrad(rho_)
                  - fvc::snGrad(pd_)
                )*mesh().magSf()
            )
        );
    }

    return tUEqn;
}


void interFluid::solvePEqn
(
    const int corr,
    const int nCorr,
    fvVectorMatrix& UEqn,
    const label pdRefCell,
    const scalar pdRefValue
)
{
    dictionary pimple = mesh().solutionDict().subDict("PIMPLE");

    volScalarField rUA = 1.0/UEqn.A();
    surfaceScalarField rUAf = fvc::interpolate(rUA);

    U() = rUA*UEqn.H();

    surfaceScalarField phiU
    (
        "phiU",
        (fvc::interpolate(U()) & mesh().Sf())
      + fvc::ddtPhiCorr(rUA, rho_, U(), phi())
    );

    adjustPhi(phiU, U(), pd_);

    phi() = phiU +
        (
            fvc::interpolate(interface_.sigmaK())*fvc::snGrad(alpha1_)
          - ghf_*fvc::snGrad(rho_)
        )*rUAf*mesh().magSf();

    const int nNonOrthCorr =
        pimple.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);

    for (int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pdEqn
        (
            fvm::laplacian(rUAf, pd_) == fvc::div(phi())
        );

        pdEqn.setReference(pdRefCell, pdRefValue);

        if (corr == nCorr - 1 && nonOrth == nNonOrthCorr)
        {
            pdEqn.solve(mesh().solutionDict().solver(pd_.name() + "Final"));
        }
        else
        {
            pdEqn.solve(mesh().solutionDict().solver(pd_.name()));
        }

        if (nonOrth == nNonOrthCorr)
        {
            phi() -= pdEqn.flux();
        }
    }

    U() += rUA*fvc::reconstruct((phi() - phiU)/rUAf);
    U().correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

interFluid::interFluid
(
    Time& runTime,
    const word& region
)
:
    fluidModel(typeName, runTime, region),
    pd_
    (
        IOobject
        (
            "pd",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
    alpha1_
    (
        IOobject
        (
            "alpha1",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
    twoPhaseProperties_(U(), phi(), "alpha1"),
    rho1_(twoPhaseProperties_.rho1()),
    rho2_(twoPhaseProperties_.rho2()),
    rho_
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT
        ),
        alpha1_*rho1_ + (scalar(1) - alpha1_)*rho2_,
        alpha1_.boundaryField().types()
    ),
    rhoPhi_
    (
        IOobject
        (
            "rho*phi",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho1_*phi()
    ),
    g_
    (
        IOobject
        (
            "g",
            runTime.constant(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    gh_("gh", g_ & mesh().C()),
    ghf_("gh", g_ & mesh().Cf()),
    interface_(alpha1_, U(), twoPhaseProperties_),
    turbulence_
    (
        incompressible::turbulenceModel::New(U(), phi(), twoPhaseProperties_)
    )
{
    // Reset p dimensions
    Info<< "Resetting the dimensions of p" << endl;
    p().dimensions().reset(dimPressure);
    p() = pd_ + rho_*gh_;

    rho_.oldTime();

    label pdRefCell = 0;
    scalar pdRefValue = 0.0;
    setRefCell(p(), fluidProperties(), pdRefCell, pdRefValue);
    //mesh().schemesDict().setFluxRequired(pd_.name());

    // Lookup pimple dict
    dictionary pimple = mesh().solutionDict().subDict("PIMPLE");

    scalar pRefValue = 0.0;

    if (pd_.needReference())
    {
        pRefValue = readScalar(pimple.lookup("pRefValue"));

        p() += dimensionedScalar
        (
            "p",
            p().dimensions(),
            pRefValue - getRefCellValue(p(), pdRefCell)
        );
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> interFluid::patchViscousForce(const label patchID) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

    tvF() =
        (
            mesh().boundary()[patchID].nf()
          & turbulence_->devReff()().boundaryField()[patchID]
        );

    vectorField n = mesh().boundary()[patchID].nf();
    tvF() -= (sqr(n) & tvF());

    return tvF;
}


tmp<scalarField> interFluid::patchPressureForce(const label patchID) const
{
    tmp<scalarField> tpF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    tpF() = pd_.boundaryField()[patchID];

    return tpF;
}


tmp<scalarField> interFluid::faceZoneMuEff
(
    const label zoneID,
    const label patchID
) const
{
    scalarField pMuEff =
       turbulence_->nuEff()().boundaryField()[patchID];

    tmp<scalarField> tMuEff
    (
        new scalarField(mesh().faceZones()[zoneID].size(), 0)
    );
    scalarField& muEff = tMuEff();

    const label patchStart =
        mesh().boundaryMesh()[patchID].start();

    forAll(pMuEff, I)
    {
        muEff[mesh().faceZones()[zoneID].whichFace(patchStart + I)] =
            pMuEff[I];
    }

    // Parallel data exchange: collect pressure field on all processors
    reduce(muEff, sumOp<scalarField>());

    return tMuEff;
}


bool interFluid::evolve()
{
    Info<< "Evolving fluid model: " << this->type() << endl;

    fvMesh& mesh = fluidModel::mesh();

    // Lookup pimple dict
    dictionary pimple = mesh.solutionDict().subDict("PIMPLE");

    // Prepare for the pressure solution
    label pdRefCell = 0;
    scalar pdRefValue = 0.0;
    setRefCell(p(), fluidProperties(), pdRefCell, pdRefValue);

    scalar pRefValue = 0.0;

    if (pd_.needReference())
    {
        pRefValue = readScalar(pimple.lookup("pRefValue"));
    }

    // Calculate CourantNo
    {
        scalar CoNum = 0.0;
        scalar meanCoNum = 0.0;
        scalar velMag = 0.0;
        CourantNo(CoNum, meanCoNum, velMag);
    }

    const int nCorr(readInt(pimple.lookup("nCorrectors")));

    const int nOuterCorr =
        readInt(fluidProperties().lookup("nOuterCorrectors"));

    // Pressure-velocity corrector
    for (int oCorr = 0; oCorr < nOuterCorr; oCorr++)
    {
        twoPhaseProperties_.correct();

        solveAlphaEqnSubCycle();

        fvVectorMatrix UEqn = solveUEqn();

        // --- PISO loop
        for (int corr = 0; corr < nCorr; corr++)
        {
            solvePEqn(corr, nCorr, UEqn, pdRefCell, pdRefValue);
        }

        fluidModel::continuityErrs();

        p() = pd_ + rho_*gh_;

        if (pd_.needReference())
        {
            p() +=
                dimensionedScalar
                (
                    "p",
                    p().dimensions(),
                    pRefValue - getRefCellValue(p(), pdRefCell)
                );
        }

        gradU() = fvc::grad(U());

        turbulence_->correct();
    }

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// ************************************************************************* //
