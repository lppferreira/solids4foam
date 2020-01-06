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

#include "buoyantBoussinesqPimpleFluid.H"
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

defineTypeNameAndDebug(buoyantBoussinesqPimpleFluid, 0);
addToRunTimeSelectionTable(physicsModel, buoyantBoussinesqPimpleFluid, fluid);
addToRunTimeSelectionTable(fluidModel, buoyantBoussinesqPimpleFluid, dictionary);


// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

void buoyantBoussinesqPimpleFluid::FourierNo()
{
    // Finds the characteristic size
    const volScalarField& cellDims = fluidModel::cellDimensions();

    FourierNum() = 0.0;
    scalar meanFourierNum = 0.0;

    FourierNum() = 
      gMax((kappaEff_ + turbulence_->nut())
           *runTime().deltaT0Value() / pow(cellDims.field(), 2));

    meanFourierNum = 
      gAverage((kappaEff_ + turbulence_->nut())
               *runTime().deltaT0Value() / pow(cellDims.field(), 2));

    Info<< "Fourier number mean: " << meanFourierNum
        << " max calculated: " << FourierNum() << endl;
}


tmp<fvVectorMatrix> buoyantBoussinesqPimpleFluid::solveUEqn()
{
    // Solve the momentum equation

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(U())
      + fvm::div(phi(), U())
      + turbulence_->divDevReff()
    );
    fvVectorMatrix& UEqn = tUEqn();

    UEqn.relax();

    if (pimple().momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                  - ghf_*fvc::snGrad(rhok_)
                  - fvc::snGrad(p_rgh_)
                )*mesh().magSf()
            )
        );
    }

    return tUEqn;
}


void buoyantBoussinesqPimpleFluid::solveTEqn()
{
    kappaEff_ =
        turbulence_->nu()/Pr_ + turbulence_->nut()/Prt_;

    fvScalarMatrix TEqn
    (
        fvm::ddt(T_)
      + fvm::div(phi(), T_)
      - fvm::laplacian(kappaEff_, T_)
    );

    TEqn.relax();

    TEqn.solve();

    rhok_ = 1.0 - beta_*(T_ - TRef_);
}


void buoyantBoussinesqPimpleFluid::solvePEqn
(
    tmp<fvVectorMatrix>& UEqn
)
{
    volScalarField rUA("rUA", 1.0/UEqn().A());
    surfaceScalarField rUAf("rUAf", fvc::interpolate(rUA));

    volVectorField HbyA("HbyA", U());
    HbyA = rUA*UEqn().H();

    surfaceScalarField phig(-rUAf*ghf_*fvc::snGrad(rhok_)*mesh().magSf());

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (fvc::interpolate(HbyA) & mesh().Sf())
      + fvc::ddtPhiCorr(rUA, U(), phi())
      + phig
    );

    adjustPhi(phiHbyA, U(), p_rgh_);

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phiHbyA, U());

    while (pimple().correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rUAf, p_rgh_) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell_, getRefCellValue(p_rgh_, pRefCell_));

        p_rghEqn.solve
        (
            mesh().solutionDict().solver(p_rgh_.select(pimple().finalInnerIter()))
        );

        if (pimple().finalNonOrthogonalIter())
        {
            // Calculate the conservative fluxes
            phi() = phiHbyA - p_rghEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p_rgh_.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U() = HbyA + rUA*fvc::reconstruct((phig - p_rghEqn.flux())/rUAf);
            U().correctBoundaryConditions();
        }
    }

    fluidModel::continuityErrs();

    p() = p_rgh_ + rhok_*gh_;

    if (p_rgh_.needReference())
    {
        p() += dimensionedScalar
        (
            "p",
            p().dimensions(),
            pRefValue_ - getRefCellValue(p(), pRefCell_)
        );
        p_rgh_ = p() - rhok_*gh_;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

buoyantBoussinesqPimpleFluid::buoyantBoussinesqPimpleFluid
(
    Time& runTime,
    const word& region
)
:
    fluidModel(typeName, runTime, region),
    p_rgh_
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
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
    laminarTransport_(U(), phi()),
    beta_(laminarTransport_.lookup("beta")),
    TRef_(laminarTransport_.lookup("TRef")),
    Pr_(laminarTransport_.lookup("Pr")),
    Prt_(laminarTransport_.lookup("Prt")),
    turbulence_
    (
        incompressible::turbulenceModel::New
        (
            U(), phi(), laminarTransport_
        )
    ),
    kappaEff_
    (
        IOobject
        (
            "kappaEff",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        turbulence_->nu()/Pr_ + turbulence_->nut()/Prt_
    ),
    rho_(laminarTransport_.lookup("rho")),
    C_(laminarTransport_.lookup("C")),
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
    hRef_
    (
        IOobject
        (
            "hRef",
            runTime.constant(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedScalar("zero", dimLength, 0)
    ),
    ghRef_(-mag(g_)*hRef_),
    gh_("gh", (g_ & mesh().C()) - ghRef_),
    ghf_("ghf", (g_ & mesh().Cf()) - ghRef_),
    rhok_
    (
        IOobject
        (
            "rhok",
            runTime.timeName(),
            mesh()
        ),
        1.0 - beta_*(T_ - TRef_)
    ),
    pRefCell_(0),
    pRefValue_(0)
{
    UisRequired();

    // Reset p dimensions
    Info<< "Resetting the dimensions of p" << endl;
    p().dimensions().reset(dimPressure/dimDensity);
    p() = p_rgh_ + rhok_*gh_;

    setRefCell(p(), pimple().dict(), pRefCell_, pRefValue_);
    mesh().schemesDict().setFluxRequired(p_rgh_.name());

    if (p_rgh_.needReference())
    {
        p() += dimensionedScalar
        (
            "p",
            p().dimensions(),
            pRefValue_ - getRefCellValue(p(), pRefCell_)
        );
        p_rgh_ = p() - rhok_*gh_;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> buoyantBoussinesqPimpleFluid::patchViscousForce
(
    const label patchID
) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

    tvF() =
        rho_.value()
       *(
            mesh().boundary()[patchID].nf()
          & (-turbulence_->devReff()().boundaryField()[patchID])
        );

    return tvF;
}


tmp<scalarField> buoyantBoussinesqPimpleFluid::patchPressureForce
(
    const label patchID
) const
{
    tmp<scalarField> tpF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    tpF() = rho_.value()*p().boundaryField()[patchID];

    return tpF;
}


tmp<scalarField> buoyantBoussinesqPimpleFluid::patchHeatFlux
(
    const label patchID
) const
{
    tmp<scalarField> thF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    thF() =
        kappaEff_.boundaryField()[patchID]
      * rho_.value()*C_.value()
      * T_.boundaryField()[patchID].snGrad();

    return thF;
}


tmp<scalarField> buoyantBoussinesqPimpleFluid::patchTemperature
(
    const label patchID
) const
{
    tmp<scalarField> tT
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    tT() = T_.boundaryField()[patchID].patchInternalField();

    return tT;
}


tmp<scalarField> buoyantBoussinesqPimpleFluid::patchKappaDelta
(
    const label patchID
) const
{
    tmp<scalarField> tKD
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    tKD() =
        kappaEff_.boundaryField()[patchID]
      * rho_.value()*C_.value()
      * mesh().boundary()[patchID].deltaCoeffs();

    return tKD;
}


bool buoyantBoussinesqPimpleFluid::evolve()
{
    Info<< "Evolving fluid model: " << this->type() << endl;

    fvMesh& mesh = fluidModel::mesh();

    bool meshChanged = false;
    if (fluidModel::fsiMeshUpdate())
    {
        // The FSI interface is in charge of calling mesh.update()
        meshChanged = fluidModel::fsiMeshUpdateChanged();
    }
    else
    {
        meshChanged = refCast<dynamicFvMesh>(mesh).update();
        reduce(meshChanged, orOp<bool>());
    }

    if (meshChanged)
    {
        const Time& runTime = fluidModel::runTime();
#       include "volContinuity.H"
    }

    // Update gh fields as the mesh may have moved
    gh_ = (g_ & mesh.C()) - ghRef_;
    ghf_ = (g_ & mesh.Cf()) - ghRef_;

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi(), U());

    // Calculate fourier number
    FourierNo();

    // Calculate CourantNo
    {
        scalar CoNum = 0.0;
        scalar meanCoNum = 0.0;
        scalar velMag = 0.0;
        CourantNo(CoNum, meanCoNum, velMag);
    }

    // Pressure-velocity corrector
    while (pimple().loop())
    {
        // Momentum equation
        tmp<fvVectorMatrix> UEqn = solveUEqn();

        // Temperature equation
        solveTEqn();

        // --- Pressure corrector loop
        while (pimple().correct())
        {
            solvePEqn(UEqn);
        }

        // Make the fluxes relative to the mesh motion
        fvc::makeRelative(phi(), U());

        UEqn.clear();

        gradU() = fvc::grad(U());

        turbulence_->correct();
    }

    Info<< "Fluid temperature min/max(T) = " << min(T_).value()
	<< ", " << max(T_).value() << " [K]" << endl;

    // Make the fluxes absolut to the mesh motion
    fvc::makeAbsolute(phi(), U());

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// ************************************************************************* //
