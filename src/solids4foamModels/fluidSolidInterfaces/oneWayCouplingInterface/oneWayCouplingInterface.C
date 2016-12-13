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

#include "oneWayCouplingInterface.H"
#include "addToRunTimeSelectionTable.H"
#include "oneWayFsiFluid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidSolidInterfaces
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(oneWayCouplingInterface, 0);
addToRunTimeSelectionTable
(
    fluidSolidInterface, oneWayCouplingInterface, dictionary
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

oneWayCouplingInterface::oneWayCouplingInterface
(
    dynamicFvMesh& fluidMesh,
    dynamicFvMesh& solidMesh
)
:
    fluidSolidInterface(typeName, fluidMesh, solidMesh)//,
//    solidZoneTraction_()
{

    // Initialize zone traction fields - JN: No need for this in of30
/*    solidZoneTraction_ =
        vectorField
        (
            solidMesh.faceZones()[solidZoneIndex()]().size(),
            vector::zero
        );*/

    if (!isA<fluidModels::oneWayFsiFluid>(fluid()))
    {
        FatalErrorIn
        (
            "oneWayCouplingInterface::oneWayCouplingInterface\n"
            "(\n"
            "    dynamicFvMesh& fluidMesh,\n"
            "    dynamicFvMesh& solidMesh\n"
            ")\n"
        )   << "The " << type() << " FSI coupling can only be used with the "
            << fluidModels::oneWayFsiFluid::typeName << " fluidModel"
            << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void oneWayCouplingInterface::evolve()
{
    fluidSolidInterface::initializeFields();

    updateInterpolator();

    fluid().evolve();

    updateTraction();

    solid().evolve();

    solid().updateTotalFields();
}


void oneWayCouplingInterface::updateTraction()
{
    Info<< "Update traction on solid patch" << endl;

    // Calculate fluid traction

    //JN: No more zones in of30

/*    const vectorField& p =
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

    vectorField fluidZoneTractionAtSolid = AMI().interpolateToTarget(-fluidZoneTraction);

    solidZoneTraction_ = fluidZoneTractionAtSolid;
*/

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

    // Solid patch total traction
    vectorField solidPatchTotalTraction =
        AMI().interpolateToTarget(-fluidPatchTotalTraction);

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
