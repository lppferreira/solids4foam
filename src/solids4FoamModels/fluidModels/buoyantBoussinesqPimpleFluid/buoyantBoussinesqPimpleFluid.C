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

tmp<fvVectorMatrix> buoyantBoussinesqPimpleFluid::solveUEqn()
{
    // Solve the momentum equation

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi(), U())
      + turbulence_->divDevReff()
    );
    fvVectorMatrix& UEqn = tUEqn();

    UEqn.relax();

    solve
    (
        UEqn
      ==
        fvc::reconstruct
        (
            (
                fvc::interpolate(rhok_)*(g_ & mesh().Sf())
              - fvc::snGrad(p())*mesh().magSf()
            )
        )
    );

    return tUEqn;
}


void buoyantBoussinesqPimpleFluid::solveTEqn()
{
    const volScalarField kappaEff
    (
        "kappaEff",
        turbulence_->nu()/Pr_ + turbulence_->nut()/Prt_
    );

    fvScalarMatrix TEqn
    (
        fvm::div(phi(), T_)
      - fvm::Sp(fvc::div(phi()), T_)
      - fvm::laplacian(kappaEff, T_)
    );

    TEqn.relax();

    TEqn.solve();

    rhok_ = 1.0 - beta_*(T_ - TRef_);
}


void buoyantBoussinesqPimpleFluid::solvePEqn
(
    pimpleControl& pimple,
    tmp<fvVectorMatrix>& UEqn,
    const label pRefCell,
    const scalar pRefValue
)
{
    volScalarField rUA("rUA", 1.0/UEqn().A());
    surfaceScalarField rUAf("(1|A(U))", fvc::interpolate(rUA));

    U() = rUA*UEqn().H();
    UEqn.clear();

    phi() = fvc::interpolate(U()) & mesh().Sf();
    adjustPhi(phi(), U(), p());

    surfaceScalarField buoyancyPhi =
        rUAf*fvc::interpolate(rhok_)*(g_ & mesh().Sf());
    phi() += buoyancyPhi;

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rUAf, p()) == fvc::div(phi())
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            // Calculate the conservative fluxes
            phi() -= pEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p().relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U() += rUA*fvc::reconstruct((buoyancyPhi - pEqn.flux())/rUAf);
            U().correctBoundaryConditions();
        }
    }

    fluidModel::continuityErrs();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

buoyantBoussinesqPimpleFluid::buoyantBoussinesqPimpleFluid
(
    Time& runTime,
    const word& region
)
:
    fluidModel(typeName, runTime, region),
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
        incompressible::RASModel::New(U(), phi(), laminarTransport_)
    ),
    rho_(laminarTransport_.lookup("rho")),
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
    betaghf_("betagh", beta_*(g_ & mesh().Cf())),
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
    FourierNo_(0.0)
{
    UisRequired();
    pisRequired();

    // Create pimple control as it is needed to set the p reference
    pimpleControl pimple(mesh());

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p(), pimple.dict(), pRefCell, pRefValue);
    mesh().schemesDict().setFluxRequired(p().name());
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const scalar& buoyantBoussinesqPimpleFluid::FourierNo() const
{
    // Finds the characteristic size
    const volScalarField& cellDims = cellDimension(runTime());

    FourierNo_ = 0.0;
    scalar meanFourierNo = 0.0;

    const volScalarField kappaEff
    (
        "kappaEff",
        turbulence_->nu()/Pr_ + turbulence_->nut()/Prt_
    );

    FourierNo_ = 
      gMax((kappaEff + turbulence_->nut())
           *runTime().deltaT0Value() / pow(cellDims.field(), 2));

    meanFourierNo = 
      gAverage((kappaEff + turbulence_->nut())
               *runTime().deltaT0Value() / pow(cellDims.field(), 2));

    Info<< "Fourier number mean: " << meanFourierNo
        << " max calculated: " << FourierNo_ << endl;

    return FourierNo_;
}


tmp<vectorField> buoyantBoussinesqPimpleFluid::patchViscousForce(const label patchID) const
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


tmp<scalarField> buoyantBoussinesqPimpleFluid::patchThermalFlux
(
    const label patchID
) const
{
    const volScalarField kappaEff
    (
        "kappaEff",
        turbulence_->nu()/Pr_ + turbulence_->nut()/Prt_
    );

    tmp<scalarField> ttF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    ttF() = kappaEff.boundaryField()[patchID]
	   * mesh().boundary()[patchID].magSf()
	   * T_.boundaryField()[patchID].snGrad();

    return ttF;
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


tmp<scalarField> buoyantBoussinesqPimpleFluid::patchKDelta
(
    const label patchID
) const
{
    const volScalarField kappaEff
    (
        "kappaEff",
        turbulence_->nu()/Pr_ + turbulence_->nut()/Prt_
    );

    tmp<scalarField> tKD
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    tKD() = kappaEff.boundaryField()[patchID]*mesh().boundary()[patchID].deltaCoeffs();

    return tKD;
}


void buoyantBoussinesqPimpleFluid::setTemperature
(
    const label patchID,
    const scalarField& temperature,
    const scalarField& nbrKDelta
)
{
    if
    (
        T_.boundaryField()[patchID].type()
     != mixedFvPatchScalarField::typeName
    )
    {
        FatalErrorIn("void interPhaseChangeFluid::setTemperature(...)")
            << "Bounary condition on " << T_.name()
                <<  " is "
                << T_.boundaryField()[patchID].type()
                << "for patch" << mesh().boundary()[patchID].name()
                << ", instead "
                << mixedFvPatchScalarField::typeName
                << abort(FatalError);
    }

    mixedFvPatchScalarField& patchT =
        refCast<mixedFvPatchScalarField>
        (
            T_.boundaryField()[patchID]
        );

    patchT.refValue() = temperature;
    patchT.refGrad() = 0.0;
    patchT.valueFraction() = nbrKDelta / (nbrKDelta + patchKDelta(patchID));
    patchT.evaluate();
}

void buoyantBoussinesqPimpleFluid::setTemperature
(
    const scalarField& faceZoneTemperature,
    const scalarField& faceZoneKDelta
)
{
    scalarField nbrpatchTemperature =
	globalPatch().globalFaceToPatch(faceZoneTemperature);

    scalarField nbrpatchKDelta =
	globalPatch().globalFaceToPatch(faceZoneKDelta);

    setTemperature(patchID, nbrpatchTemperature, nbrpatchKDelta);
}


bool buoyantBoussinesqPimpleFluid::evolve()
{
    Info<< "Evolving fluid model: " << this->type() << endl;

    fvMesh& mesh = fluidModel::mesh();

    // Create pimple control
    pimpleControl pimple(mesh);

    // Prepare for the pressure solution
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p(), pimple.dict(), pRefCell, pRefValue);

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
    while (pimple.loop())
    {
        // Momentum equation
        tmp<fvVectorMatrix> UEqn = solveUEqn();

        // Temperature equation
        solveTEqn();

        // Pressure equation
        solvePEqn(pimple, UEqn, pRefCell, pRefValue);

        gradU() = fvc::grad(U());

        turbulence_->correct();
    }

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// ************************************************************************* //
