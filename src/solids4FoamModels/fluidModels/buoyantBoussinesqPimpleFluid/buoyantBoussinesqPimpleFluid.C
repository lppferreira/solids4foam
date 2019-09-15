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
#include "mixedFvPatchFields.H"

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
                    fvc::interpolate(rhok_)*(g_ & mesh().Sf())
                  - fvc::snGrad(p())*mesh().magSf()
                )
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
        fvm::ddt(T())
      + fvm::div(phi(), T())
      - fvm::laplacian(kappaEff_, T())
    );

    TEqn.relax();

    TEqn.solve();

    rhok_ = 1.0 - beta_*(T() - TRef_);
}


void buoyantBoussinesqPimpleFluid::solvePEqn
(
    tmp<fvVectorMatrix>& UEqn
)
{
    volScalarField rUA("rUA", 1.0/UEqn().A());
    surfaceScalarField rUAf("(1|A(U))", fvc::interpolate(rUA));

    U() = rUA*UEqn().H();

    phi() = fvc::interpolate(U()) & mesh().Sf();
    adjustPhi(phi(), U(), p());

    surfaceScalarField buoyancyPhi =
        rUAf*fvc::interpolate(rhok_)*(g_ & mesh().Sf());
    phi() += buoyancyPhi;

    while (pimple().correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rUAf, p()) == fvc::div(phi())
        );

        pEqn.setReference(pRefCell_, pRefValue_);

        pEqn.solve();

        if (pimple().finalNonOrthogonalIter())
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
    betaghf_("betagh", beta_*(g_ & mesh().Cf())),
    rhok_
    (
        IOobject
        (
            "rhok",
            runTime.timeName(),
            mesh()
        ),
        1.0 - beta_*(T() - TRef_)
    ),
    pRefCell_(0),
    pRefValue_(0),
    FourierNo_(0)
{
    UisRequired();
    pisRequired();
    TisRequired();

    setRefCell(p(), pimple().dict(), pRefCell_, pRefValue_);
    mesh().schemesDict().setFluxRequired(p().name());
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar& buoyantBoussinesqPimpleFluid::FourierNo()
{
    // Finds the characteristic size
    const volScalarField& cellDims = cellDimension(runTime());

    FourierNo_ = 0.0;
    scalar meanFourierNo = 0.0;

    FourierNo_ = 
      gMax((kappaEff_ + turbulence_->nut())
           *runTime().deltaT0Value() / pow(cellDims.field(), 2));

    meanFourierNo = 
      gAverage((kappaEff_ + turbulence_->nut())
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
    tmp<scalarField> ttF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    ttF() = fvc::interpolate(kappaEff_)().boundaryField()[patchID]
          * rho_.value() * C_.value()
          * mesh().boundary()[patchID].magSf()
          * T().boundaryField()[patchID].snGrad();

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

    tT() = T().boundaryField()[patchID].patchInternalField();

    return tT;
}


tmp<scalarField> buoyantBoussinesqPimpleFluid::patchKDelta
(
    const label patchID
) const
{
    tmp<scalarField> tKD
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    tKD() = fvc::interpolate(kappaEff_)().boundaryField()[patchID]
          * rho_.value() * C_.value()
          * mesh().boundary()[patchID].deltaCoeffs();

    return tKD;
}


void buoyantBoussinesqPimpleFluid::setTemperature
(
    const label patchID,
    const scalarField& faceZoneTemperature,
    const scalarField& faceZoneKDelta
)
{
    if
    (
        T().boundaryField()[patchID].type()
     != mixedFvPatchScalarField::typeName
    )
    {
        FatalErrorIn("void buoyantBoussinesqPimpleFluid::setTemperature(...)")
            << "Bounary condition on " << T().name()
                <<  " is "
                << T().boundaryField()[patchID].type()
                << "for patch" << mesh().boundary()[patchID].name()
                << ", instead of "
                << mixedFvPatchScalarField::typeName
                << abort(FatalError);
    }

    scalarField nbrPatchTemperature =
	globalPatch().globalFaceToPatch(faceZoneTemperature);

    scalarField nbrPatchKDelta =
	globalPatch().globalFaceToPatch(faceZoneKDelta);

    mixedFvPatchScalarField& patchT =
        refCast<mixedFvPatchScalarField>
        (
            T().boundaryField()[patchID]
        );

    patchT.refValue() = nbrPatchTemperature;
    patchT.refGrad() = 0.0;
    patchT.valueFraction() = nbrPatchKDelta / (nbrPatchKDelta + patchKDelta(patchID));
    patchT.evaluate();
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

        // Pressure equation
        while (pimple().correct())
        {
            solvePEqn(UEqn);
        }

        UEqn.clear();

        // Make the fluxes relative to the mesh motion
        fvc::makeRelative(phi(), U());

        gradU() = fvc::grad(U());

        turbulence_->correct();
    }

    // Make the fluxes absolut to the mesh motion
    fvc::makeAbsolute(phi(), U());

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// ************************************************************************* //
