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

#include "fluidSolidInterface.H"
#include "volFields.H"
#include "polyPatchID.H"
#include "primitivePatchInterpolation.H"
#include "twoDPointCorrector.H"
#include "fixedValuePointPatchFields.H"
#include "ZoneIDs.H"
#include "SubField.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fluidSolidInterface, 0);
    defineRunTimeSelectionTable(fluidSolidInterface, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::fluidSolidInterface::updateCoupled()
{
    if (couplingStartTime_ > SMALL && !coupled_)
    {
        if (runTime().value() > (couplingStartTime_ - SMALL))
        {
            InfoIn("fluidSolidInterface::updateCoupled()")
                << "Enabling fluid-solid coupling" << endl;

            // Enable coupling
            coupled_ = true;

            return true;
        }
    }

    return false;
}


void Foam::fluidSolidInterface::calcAMIInterpolator() const
{
    // Create AMI interpolation
    if (!AMIPtr_.empty())
    {
        FatalErrorIn
        (
            "void fluidSolidInterface::calcAMIInterpolator() const"
        )   << "pointer already exists"
            << abort(FatalError);
    }

    AMIPtr_.set
    (
        new AMIPatchToPatchInterpolation
        (
            fluidMesh().boundaryMesh()[fluidPatchIndex_],
            solidMesh().boundaryMesh()[solidPatchIndex_],
            faceAreaIntersect::tmMesh,
            true,
            AMIPatchToPatchInterpolation::imFaceAreaWeight,
            -1,
            false // flip
        )
    );

    Info<< "Checking fluid-to-solid face interpolator (AMI)" << endl;
    {
        const vectorField& fluidPatchFaceCentres =
            fluidMesh().boundaryMesh()[fluidPatchIndex_].faceCentres();

        const vectorField solidPatchFaceCentres =
            AMI().interpolateToTarget(fluidPatchFaceCentres);

        const scalar maxDist =
            gMax
            (
                mag
                (
                    solidPatchFaceCentres
                  - solidMesh().boundaryMesh()[solidPatchIndex_].faceCentres()
                )
            );

        Info<< "    Fluid-to-solid face interpolation error: " << maxDist
            << endl;
    }
}


void Foam::fluidSolidInterface::calcAccumulatedFluidInterfaceDisplacement() const
{
    // Read accumulated displacement
    if (accumulatedFluidInterfaceDisplacementPtr_)
    {
        FatalErrorIn
        (
            "void fluidSolidInterface::"
            "calcAccumulatedFluidInterfaceDisplacement() const"
        )   << "Accumulated displacement field already exists"
            << abort(FatalError);
    }

    // Accumulated fluid interface displacement
    IOobject accumulatedFluidInterfaceDisplacementHeader
    (
        "accumulatedFluidInterfaceDisplacement",
        fluid().runTime().timeName(),
        fluidMesh(),
        IOobject::MUST_READ
    );

    if (accumulatedFluidInterfaceDisplacementHeader.headerOk())
    {
        Pout << "Reading accumulated fluid interface displacement" << endl;

        accumulatedFluidInterfaceDisplacementPtr_ =
            new vectorIOField
            (
                IOobject
                (
                    "accumulatedFluidInterfaceDisplacement",
                    fluid().runTime().timeName(),
                    fluidMesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )
            );
    }
    else
    {
        Pout<< "Creating accumulated fluid interface displacement" << endl;

        accumulatedFluidInterfaceDisplacementPtr_ =
            new vectorIOField
            (
                IOobject
                (
                    "accumulatedFluidInterfaceDisplacement",
                    fluid().runTime().timeName(),
                    fluidMesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                vectorField
                (
                    fluidMesh().boundaryMesh()[fluidPatchIndex()].nPoints(),
                    vector::zero
                )
            );
    }
}


Foam::vectorIOField&
Foam::fluidSolidInterface::accumulatedFluidInterfaceDisplacement()
{
    if (!accumulatedFluidInterfaceDisplacementPtr_)
    {
        calcAccumulatedFluidInterfaceDisplacement();
    }

    return *accumulatedFluidInterfaceDisplacementPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidSolidInterface::fluidSolidInterface
(
    const word& type,
    dynamicFvMesh& fluidMesh,
    dynamicFvMesh& solidMesh
)
:
    IOdictionary
    (
        IOobject
        (
            "fsiProperties",
            fluidMesh.time().constant(),
            fluidMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    fsiProperties_(subDict(type + "Coeffs")),
    fluidMesh_(fluidMesh),
    fluid_(fluidModel::New(fluidMesh_)),
    solidMesh_(solidMesh),
    solid_(solidModel::New(solidMesh_)),
    solidPatchIndex_(-1),
    fluidPatchIndex_(-1),
    AMIPtr_(NULL),
    outerCorrTolerance_
    (
        fsiProperties_.lookupOrDefault<scalar>("outerCorrTolerance", 1e-06)
    ),
    nOuterCorr_
    (
        fsiProperties_.lookupOrDefault<int>("nOuterCorr", 30)
    ),
    coupled_
    (
        fsiProperties_.lookupOrDefault<Switch>("coupled", true)
    ),
    couplingStartTime_
    (
        fsiProperties_.lookupOrDefault<scalar>("couplingStartTime", -1.0)
    ),
    interfaceDeformationLimit_
    (
        fsiProperties_.lookupOrDefault<scalar>("interfaceDeformationLimit", 0.0)
    ),
    fluidPatchPointsDispl_(),
    fluidPatchPointsDisplRef_(),
    fluidPatchPointsDisplPrev_(),
    solidPatchPointsDispl_(),
    solidPatchPointsDisplRef_(),
    residual_(),
    residualPrev_(),
    maxResidualNorm_(0),
    outerCorr_(0),
    interpolatorUpdateFrequency_
    (
        fsiProperties_.lookupOrDefault<int>("interpolatorUpdateFrequency", 0)
    ),
    accumulatedFluidInterfaceDisplacementPtr_(NULL)
{
    // Check if couplingStartTime is specified
    if (couplingStartTime_ > SMALL)
    {
        if (coupled_)
        {
            WarningIn(type + "::fsiProperties(...)")
                << "When using the coupilngStartTime option, the coupled "
                << "option should be set to off: resetting coupled to off"
                 << endl;

            coupled_ = false;
        }
    }

    // Solid patch index

    const word solidPatchName(fsiProperties_.lookup("solidPatch"));

    const polyPatchID solidPatch
    (
        solidPatchName,
        solidMesh.boundaryMesh()
    );

    if (!solidPatch.active())
    {
        FatalErrorIn("fluidSolidInterface::fluidSolidInterface(...)")
            << "Solid patch name " << solidPatchName << " not found."
            << abort(FatalError);
    }

    Info<< "Reading solidPatch: " << solidPatchName << endl;
    solidPatchIndex_ = solidPatch.index();

    // Fluid patch index

    const word fluidPatchName(fsiProperties_.lookup("fluidPatch"));

    const polyPatchID fluidPatch
    (
        fluidPatchName,
        fluidMesh.boundaryMesh()
    );

    if (!fluidPatch.active())
    {
        FatalErrorIn("fluidSolidInterface::fluidSolidInterface(...)")
            << "Fluid patch name " << fluidPatchName << " not found."
            << abort(FatalError);
    }

    Info<< "Reading fluidPatch: " << fluidPatchName << endl;
    fluidPatchIndex_ = fluidPatch.index();

    // Philip: why do we initialise the residual here, but we initialise other
    // interface fields in the initialiseFields function?
    // Initialize residual
    residual_ =
        vectorField
        (
            fluidMesh.boundaryMesh()[fluidPatchIndex_].nPoints(),
            vector::zero
        );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidSolidInterface::~fluidSolidInterface()
{
    AMIPtr_.clear();
    deleteDemandDrivenData(accumulatedFluidInterfaceDisplacementPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


const Foam::AMIPatchToPatchInterpolation& Foam::fluidSolidInterface::AMI() const
{
    if (!AMIPtr_.valid())
    {
        calcAMIInterpolator();
    }

    return AMIPtr_;
}


void Foam::fluidSolidInterface::initializeFields()
{
    fluidPatchPointsDispl_ =
        vectorField
        (
            fluidMesh().boundaryMesh()[fluidPatchIndex_].nPoints(),
            vector::zero
        );

    fluidPatchPointsDisplRef_ =
        vectorField
        (
            fluidMesh().boundaryMesh()[fluidPatchIndex_].nPoints(),
            vector::zero
        );

    fluidPatchPointsDisplPrev_ =
        vectorField
        (
            fluidMesh().boundaryMesh()[fluidPatchIndex_].nPoints(),
            vector::zero
        );

    solidPatchPointsDispl_ =
        vectorField
        (
            fluidMesh().boundaryMesh()[fluidPatchIndex_].nPoints(),
            vector::zero
        );

    solidPatchPointsDisplRef_ =
        vectorField
        (
            fluidMesh().boundaryMesh()[fluidPatchIndex_].nPoints(),
            vector::zero
        );

    residualPrev_ = residual_;

    residual_ =
        vectorField
        (
            fluidMesh().boundaryMesh()[fluidPatchIndex_].nPoints(),
            vector::zero
        );

    maxResidualNorm_ = 0;

    outerCorr_ = 0;
}


void Foam::fluidSolidInterface::updateInterpolator()
{
    if (interpolatorUpdateFrequency_ != 0)
    {
        if (((runTime().timeIndex() - 1) % interpolatorUpdateFrequency_) == 0)
        {
            AMIPtr_.clear();
            AMI();
        }
    }
    else
    {
        if ((runTime().timeIndex() - 1) == 0)
        {
            AMIPtr_.clear();
            AMI();
        }
    }
}


void Foam::fluidSolidInterface::moveFluidMesh()
{
    // Get fluid patch displacements
    const vectorField& fluidPatchPointsDispl = this->fluidPatchPointsDispl();
    const vectorField& fluidPatchPointsDisplPrev =
        this->fluidPatchPointsDisplPrev();

    // Move fluid mesh
    const vectorField& n =
        fluidMesh().boundaryMesh()[fluidPatchIndex()].pointNormals();

    // Philip: is it OK to use primitivePatchInterpolation for faceToPoint in
    // parallel? It will not be correct on shared points!
    // We could use pointVectorField.correctBoundaryConditions() to sync the
    // shared points...
    primitivePatchInterpolation patchInterpolator
    (
        fluidMesh().boundaryMesh()[fluidPatchIndex()]
    );

    const scalarField pointDeltaCoeffs =
        patchInterpolator.faceToPointInterpolate
        (
            fluidMesh().boundary()[fluidPatchIndex()].deltaCoeffs()
        );

    const scalar delta =
        gMax
        (
            mag
            (
                n
              & (
                    accumulatedFluidInterfaceDisplacement()
                  + fluidPatchPointsDispl
                  - fluidPatchPointsDisplPrev
                )
            )*pointDeltaCoeffs
        );

    Info<< "Maximal accumulated displacement of interface points: "
        << delta << endl;

    if (delta < interfaceDeformationLimit())
    {
        // Move only interface points
        pointField newPoints = fluidMesh().points();

        const labelList& meshPoints =
            fluidMesh().boundaryMesh()[fluidPatchIndex()].meshPoints();

        forAll(fluidPatchPointsDispl, pointI)
        {
            newPoints[meshPoints[pointI]] +=
                fluidPatchPointsDispl[pointI]
              - fluidPatchPointsDisplPrev[pointI];
        }

        twoDPointCorrector twoDCorrector(fluidMesh());

        twoDCorrector.correctPoints(newPoints);

        fluidMesh_.movePoints(newPoints);

        // Accumulate interface points displacement
        accumulatedFluidInterfaceDisplacement() +=
            fluidPatchPointsDispl
          - fluidPatchPointsDisplPrev;
    }
    else
    {
        // Philip comment:
        // We should consider a nicer way to allow use of different fluid mesh
        // motion solvers, instead of hard-coding if-else-if loops

        // Move whole fluid mesh
        pointField newPoints = fluidMesh().points();

        const labelList& meshPoints =
            fluidMesh().boundaryMesh()[fluidPatchIndex()].meshPoints();

        forAll (accumulatedFluidInterfaceDisplacement(), pointI)
        {
            newPoints[meshPoints[pointI]] -=
                accumulatedFluidInterfaceDisplacement()[pointI];
        }

        twoDPointCorrector twoDCorrector(fluidMesh());

        twoDCorrector.correctPoints(newPoints);

        fluidMesh_.movePoints(newPoints);

        accumulatedFluidInterfaceDisplacement() +=
            fluidPatchPointsDispl
          - fluidPatchPointsDisplPrev;

        // Check mesh motion solver type
        // bool feMotionSolver =
        //     fluidMesh().objectRegistry::foundObject<tetPointVectorField>
        //     (
        //         "motionU"
        //     );

        bool fvMotionSolver =
            fluidMesh().objectRegistry::foundObject<pointVectorField>
            (
                "pointMotionU"
            );

        // bool rbfMotionSolver =
        //     fluidMesh().objectRegistry::foundObject<RBFMotionSolver>
        //     (
        //         "dynamicMeshDict"
        //     );

//         if (rbfMotionSolver)
//         {
//             // Grab RBF motion solver
//             RBFMotionSolver& ms =
//                 const_cast<RBFMotionSolver&>
//                 (
//                     fluidMesh().objectRegistry::lookupObject<RBFMotionSolver>
//                     (
//                         "dynamicMeshDict"
//                     )
//                 );

//             Info<< "RBF mesh motion" << endl;

//             const labelList& movingMeshPoints = ms.movingIDs();

//             vectorField motion(movingMeshPoints.size(), vector::zero);

//             vectorField fluidPatchDisplacement =
//                 accumulatedFluidInterfaceDisplacement();
// //                /fluid().runTime().deltaT().value();

//             const labelList& meshPoints =
//                 fluidMesh().boundaryMesh()[fluidPatchIndex()].meshPoints();

//             forAll(meshPoints, pointI)
//             {
//                 label curMovingPoint =
//                     findIndex(movingMeshPoints, meshPoints[pointI]);

//                 if (curMovingPoint != -1)
//                 {
//                     motion[curMovingPoint] = fluidPatchDisplacement[pointI];
//                 }
//             }

//             ms.setMotion(motion);
//         }
//         else if (feMotionSolver)
//        if (feMotionSolver)
//        {
//            tetPointVectorField& motionU =
//                const_cast<tetPointVectorField&>
//                (
//                    fluidMesh().objectRegistry::
//                    lookupObject<tetPointVectorField>
//                    (
//                        "motionU"
//                    )
//                );

//            fixedValueTetPolyPatchVectorField& motionUFluidPatch =
//                refCast<fixedValueTetPolyPatchVectorField>
//                (
//                    motionU.boundaryField()[fluidPatchIndex()]
//                );

//            tetPolyPatchInterpolation tppi
//            (
//                refCast<const faceTetPolyPatch>(motionUFluidPatch.patch())
//            );

//            motionUFluidPatch ==
//                tppi.pointToPointInterpolate
//                (
//                    accumulatedFluidInterfaceDisplacement()
//                   /fluid().runTime().deltaT().value()
//                );
//        }
//        else if (fvMotionSolver)
        if (fvMotionSolver)
        {
            Info<< "Using fvMotionSolver for movement of fluid mesh..."
                << endl;

            pointVectorField& motionU =
                const_cast<pointVectorField&>
                (
                    fluidMesh().objectRegistry::
                    lookupObject<pointVectorField>
                    (
                        "pointMotionU"
                    )
                );

            fixedValuePointPatchVectorField& motionUFluidPatch =
                refCast<fixedValuePointPatchVectorField>
                (
                    motionU.boundaryField()[fluidPatchIndex()]
                );

            motionUFluidPatch ==
                accumulatedFluidInterfaceDisplacement()
               /fluid().runTime().deltaT().value();

            InfoIn("moveFluidMesh")
                << "motionUFluidPatch: "
                << max(mag(motionUFluidPatch)) << endl;
        }
        else
        {
            FatalErrorIn("fluidSolidInterface::moveFluidMesh()")
                << "Problem with fluid mesh motion solver selection"
                << abort(FatalError);
        }

        fluidMesh_.update();

        accumulatedFluidInterfaceDisplacement() =
            vectorField
            (
                accumulatedFluidInterfaceDisplacement().size(),
                vector::zero
            );
    }
}


void Foam::fluidSolidInterface::updateForce()
{
    Info<< "Setting traction on solid patch" << endl;

    // Philip: there are no face zones on of30 branch, so we instead directly
    // set patch tractions

    const vectorField fluidPatchTraction =
        fluid().patchViscousForce(fluidPatchIndex());

    const scalarField fluidPatchPressure =
        fluid().patchPressureForce(fluidPatchIndex());

    // Fluid patch face normals
    const vectorField n = fluidMesh().boundary()[fluidPatchIndex_].nf();

    // Fluid patch total traction
    const vectorField fluidPatchTotalTraction =
        fluidPatchTraction - fluidPatchPressure*n;

    // Philip: What is this used for?
    // Solid patch traction
    // vectorField solidPatchTraction =
    //     AMI().interpolateToTarget(-fluidPatchTraction);

    // Solid patch total traction
    vectorField solidPatchTotalTraction =
        AMI().interpolateToTarget(-fluidPatchTotalTraction);

    // const scalarField solidZoneMuEff =
    //     AMI().interpolateToTarget
    //     (
    //         fluid().faceZoneMuEff(fluidZoneIndex(), fluidPatchIndex())
    //     );

    // const tensorField solidZoneSurfaceGradientOfVelocity =
    //     solid().faceZoneSurfaceGradientOfVelocity
    //     (
    //         solidZoneIndex(),
    //         solidPatchIndex()
    //     );

    // Philip: What is this used for?
    // Add part of the viscous force present only
    // at the deforming and moving walls
    // solidZoneTraction +=
    //     solidZoneMuEff
    //    *(
    //        -2*tr(solidZoneSurfaceGradientOfVelocity)*solidZoneNormal
    //       + (solidZoneSurfaceGradientOfVelocity&solidZoneNormal)
    //     );

    // Philip: What is this used for?
    // const vectorField movingTraction =
    //     solidZoneMuEff
    //    *(
    //        -2*tr(solidZoneSurfaceGradientOfVelocity)*solidZoneNormal
    //       + (solidZoneSurfaceGradientOfVelocity&solidZoneNormal)
    //     );

    // Philip: What is this used for?
    //solidZonePressure_ = AMI().interpolateToTarget(fluidZonePressure);

    if (!coupled_)
    {
        updateCoupled();
    }

    if (coupled())
    {
        // No face zones on of30 branch, so we instead directly set patch
        // tractions
        solid().setTraction
        (
            solidPatchIndex(),
            solidPatchTotalTraction
        );
    }

    // Total force at the fluid side of the interface
    {
        const scalarField& magSf =
            fluidMesh().boundary()[fluidPatchIndex_].magSf();

        const vector totalTractionForce = gSum(fluidPatchTotalTraction*magSf);

        Info<< "Total force (fluid) = " << totalTractionForce << endl;
    }

    // Total force at the solid side of the interface
    {
        const scalarField& magSf =
            solidMesh().boundary()[solidPatchIndex_].magSf();

        const vector totalTractionForce = gSum(solidPatchTotalTraction*magSf);

        Info<< "Total force (solid) = " << totalTractionForce << endl;
    }
}


Foam::scalar Foam::fluidSolidInterface::updateResidual()
{
    // 1. use AMI to interpolate the solid patch face D field to the fluid
    //    patch faces D
    // 2. use primitivePatchInterpolation to interpolate the fluid patch faces D
    //    to the fluid patch points D
    // 3. assemble a pointD field and call correctBoundaryConditions to sync
    //    shared points

    // Get increment of displacement at solid interface faces
    const vectorField solidPatchDisplAtSolid =
        solid().patchDisplacementIncrement(solidPatchIndex());

    // Interpolate solid displacements to the fluid interface
    const vectorField solidPatchDispl =
        AMI().interpolateToSource(solidPatchDisplAtSolid);

    // Create patch interpolator
    // Note: the global shared points are not treated correctly in parallel
    // so we will sync them below
    primitivePatchInterpolation patchInterp
    (
        fluidMesh().boundaryMesh()[fluidPatchIndex_]
    );

    // Interpolate fluid face values to points
    solidPatchPointsDispl() = patchInterp.faceToPointInterpolate(solidPatchDispl);

    // Create pointMesh
    pointMesh pMesh(fluidMesh());

    // Create pointDD field to sync global points
    pointVectorField pointDD
    (
        IOobject
        (
            "fluidPointDD",
            fluidMesh().time().timeName(),
            fluidMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh,
        dimensionedVector("zero", dimLength, vector::zero),
        "calculated"
    );

    // Populate pointDD with the interface values
    const labelList& meshPoints =
        fluidMesh().boundaryMesh()[fluidPatchIndex_].meshPoints();

    vectorField& pointDDI = pointDD.internalField();
    vectorField& solidPatchPointsDispl = this->solidPatchPointsDispl();

    forAll(meshPoints, pI)
    {
        pointDDI[meshPoints[pI]] = solidPatchPointsDispl[pI];
    }

    // Sync the pointDD global shared points
    pointDD.correctBoundaryConditions();

    // Update the point patch displacements
    forAll(meshPoints, pI)
    {
        solidPatchPointsDispl[pI] = pointDDI[meshPoints[pI]];
    }

    // Update residual

    residualPrev() = residual();

    residual() = solidPatchPointsDispl - fluidPatchPointsDispl();

    scalar residualNorm = ::sqrt(gSum(magSqr(residual())));

    if (residualNorm > maxResidualNorm_)
    {
        maxResidualNorm_ = residualNorm;
    }

    residualNorm /= maxResidualNorm_ + SMALL;

    Info<< "Current FSI relative residual norm: " << residualNorm << endl;

    return residualNorm;
}


// ************************************************************************* //
