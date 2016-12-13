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

#include "weakCouplingInterface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidSolidInterfaces
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(weakCouplingInterface, 0);
addToRunTimeSelectionTable
(
    fluidSolidInterface, weakCouplingInterface, dictionary
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

weakCouplingInterface::weakCouplingInterface
(
    dynamicFvMesh& fluidMesh,
    dynamicFvMesh& solidMesh
)
:
    fluidSolidInterface(typeName, fluidMesh, solidMesh),
//    solidZoneTraction_(),
//    solidZoneTractionPrev_(),
    predictedSolidZoneTraction_(),//JN: Do we need this?
    relaxationFactor_
    (
        fsiProperties().lookupOrDefault<scalar>("relaxationFactor", 0.01)
    )
{
    // Initialize zone traction fields - no need for this in of30
/*    solidZoneTraction_ =
        vectorField
        (
            solidMesh.faceZones()[solidZoneIndex()]().size(),
            vector::zero
        );

    solidZoneTractionPrev_ =
        vectorField
        (
            solidMesh.faceZones()[solidZoneIndex()]().size(),
            vector::zero
        );*/
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void weakCouplingInterface::evolve()
{
    initializeFields();

    updateInterpolator();

    solid().evolve();

    updateWeakDisplacement();

    moveFluidMesh();

    fluid().evolve();

    updateWeakTraction();

    solid().updateTotalFields();
}


void weakCouplingInterface::initializeFields()
{
    predictedSolidZoneTraction_ =
        vectorField
        (
            fluidMesh().boundaryMesh()[fluidPatchIndex()].nPoints(),
            vector::zero
        );

    fluidSolidInterface::initializeFields();
}


void weakCouplingInterface::updateWeakDisplacement()
{
/*    vectorField solidZonePointsDisplAtSolid =
        solid().faceZonePointDisplacementIncrement(solidZoneIndex());

    solidZonePointsDispl() = AMI().pointInterpolateToSource(solidZonePointsDisplAtSolid);

    residualPrev() = residual();

    residual() = solidZonePointsDispl() - fluidZonePointsDispl();

    fluidZonePointsDisplPrev() = fluidZonePointsDispl();

    fluidZonePointsDispl() += residual();*/

    fluidPatchPointsDisplPrev() = fluidPatchPointsDispl();

    fluidPatchPointsDispl() += relaxationFactor_*residual();

    // Make sure that displacement on all processors is equal to one
    // calculated on master processor
/*    if (Pstream::parRun())
    {
        if(!Pstream::master())
        {
            fluidZonePointsDispl() *= 0.0;
        }

        //- pass to all procs
        reduce(fluidZonePointsDispl(), sumOp<vectorField>());

        label globalFluidZoneIndex =
            findIndex(fluid().globalFaceZones(), fluidZoneIndex());

        if (globalFluidZoneIndex == -1)
        {
            FatalErrorIn
            (
                "void weakCouplingInterface::updateWeakDisplacement()"
            )   << "global zone point map is not available"
                << abort(FatalError);
        }

        const labelList& map =
            fluid().globalToLocalFaceZonePointMap()[globalFluidZoneIndex];

        if (!Pstream::master())
        {
            vectorField fluidZonePointsDisplGlobal =
                fluidZonePointsDispl();

            forAll(fluidZonePointsDisplGlobal, globalPointI)
            {
                label localPoint = map[globalPointI];

                fluidZonePointsDispl()[localPoint] =
                    fluidZonePointsDisplGlobal[globalPointI];
            }
        }
    }*/
}


void weakCouplingInterface::updateWeakTraction()
{
    Info<< "Update weak traction on solid patch" << endl;

    //JN: No more zones in of30

/*    solidZoneTractionPrev_ = solidZoneTraction_;

    // Calc fluid traction

    const vectorField& p =
        fluidMesh().faceZones()[fluidZoneIndex()]().localPoints();
    const faceList& f =
        fluidMesh().faceZones()[fluidZoneIndex()]().localFaces();

    vectorField n(f.size(), vector::zero);
    forAll(n, faceI)
    {
        n[faceI] = f[faceI].normal(p);
        n[faceI] /= mag(n[faceI]);
    }

    vectorField fluidZoneTraction =
        fluid().faceZoneViscousForce
        (
            fluidZoneIndex(),
            fluidPatchIndex()
        )
      - fluid().faceZonePressureForce(fluidZoneIndex(), fluidPatchIndex())*n;

        vectorField fluidZoneTractionAtSolid =
            AMI().interpolateToTarget(-fluidZoneTraction);

    solidZoneTraction_ =
        relaxationFactor_*fluidZoneTractionAtSolid
      + (1.0 - relaxationFactor_)*predictedSolidZoneTraction_;*/

    //JN: instead of zones use the patches

    const vectorField fluidPatchTraction =
        fluid().patchViscousForce(fluidPatchIndex());

    const scalarField fluidPatchPressure =
        fluid().patchPressureForce(fluidPatchIndex());

    // Fluid patch face normals
    const vectorField n = fluidMesh().boundary()[fluidPatchIndex()].nf();

    // Fluid patch total traction
    const vectorField fluidPatchTotalTraction =
        fluidPatchTraction - fluidPatchPressure*n;

    //JN: Where do we update predictedSolidZoneTraction?
    //If it is zero,  we could delete it.
    // Solid patch total traction
    vectorField solidPatchTotalTraction =
        relaxationFactor_*AMI().interpolateToTarget(-fluidPatchTotalTraction);
//	+ (1.0 - relaxationFactor_)*predictedSolidZoneTraction_;

    if (coupled())
    {
        solid().setTraction
        (
            solidPatchIndex(),
            solidPatchTotalTraction
        );
    }

    // Total force at the fluid side of the interface
    {
        const scalarField& magSf =
            fluidMesh().boundary()[fluidPatchIndex()].magSf();

        const vector totalTractionForce = gSum(fluidPatchTotalTraction*magSf);

        Info<< "Total force (fluid) = " << totalTractionForce << endl;
    }

    // Total force at the solid side of the interface
    {
        const scalarField& magSf =
            solidMesh().boundary()[solidPatchIndex()].magSf();

        const vector totalTractionForce = gSum(solidPatchTotalTraction*magSf);

        Info<< "Total force (solid) = " << totalTractionForce << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolidInterfaces

} // End namespace Foam

// ************************************************************************* //
