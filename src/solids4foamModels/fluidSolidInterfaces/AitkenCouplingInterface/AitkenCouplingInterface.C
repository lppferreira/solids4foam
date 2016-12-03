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

#include "AitkenCouplingInterface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidSolidInterfaces
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(AitkenCouplingInterface, 0);
addToRunTimeSelectionTable
(
    fluidSolidInterface, AitkenCouplingInterface, dictionary
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AitkenCouplingInterface::AitkenCouplingInterface
(
    dynamicFvMesh& fluidMesh,
    dynamicFvMesh& solidMesh
)
:
    fluidSolidInterface(typeName, fluidMesh, solidMesh),
    relaxationFactor_
    (
        fsiProperties().lookupOrDefault<scalar>("relaxationFactor", 0.01)
    ),
    aitkenRelaxationFactor_(relaxationFactor_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void AitkenCouplingInterface::evolve()
{
    Info<< "initializeFields" << endl;
    initializeFields();

    Info<< "updateInterpolator" << endl;
    updateInterpolator();

    scalar residualNorm = 0;

    do
    {
        outerCorr()++;

        // Transfer the displacement from the solid to the fluid
        Info<< "updateDisplacement" << endl;
        updateDisplacement();

        // Move the fluid mesh
        Info<< "moveFluidMesh" << endl;
        moveFluidMesh();

        // Solve fluid
        Info<< "fluid.evolve" << endl;
        fluid().evolve();

        // Transfer the force from the fluid to the solid
        Info<< "updateForce" << endl;
        updateForce();

        // Solve solid
        Info<< "solid.evolve" << endl;
        solid().evolve();

        // Calculate the FSI residual
        Info<< "updateResidual" << endl;
        residualNorm = updateResidual();
        Info<< "end of loop" << endl;
    }
    while (residualNorm > outerCorrTolerance() && outerCorr() < nOuterCorr());

    solid().updateTotalFields();
}


void AitkenCouplingInterface::updateDisplacement()
{
    Info<< nl << "Time = " << fluid().runTime().timeName()
        << ", iteration: " << outerCorr() << endl;

    if (outerCorr() < 3)
    {
        Info<< "Current fsi under-relaxation factor: "
            << relaxationFactor_ << endl;

        fluidPatchPointsDisplPrev() = fluidPatchPointsDispl();

        fluidPatchPointsDispl() += relaxationFactor_*residual();
    }
    else
    {
        aitkenRelaxationFactor_ =
            -aitkenRelaxationFactor_
           *(
                sum
                (
                    residualPrev()
                  & (residual() - residualPrev())
                )
               /(
                    sum
                    (
                        (residual() - residualPrev())
                      & (residual() - residualPrev())
                    )
                )
            );

        if (Pstream::parRun())
        {
            if (!Pstream::master())
            {
                aitkenRelaxationFactor_ = 0.0;
            }

            // Pass to all procs
            reduce(aitkenRelaxationFactor_, sumOp<scalar>());
        }

        aitkenRelaxationFactor_ = mag(aitkenRelaxationFactor_);

        if (aitkenRelaxationFactor_ > 1)
        {
            aitkenRelaxationFactor_ = relaxationFactor_;
        }

        Info<< "Current fsi under-relaxation factor (Aitken): "
            << aitkenRelaxationFactor_ << endl;

        fluidPatchPointsDisplPrev() = fluidPatchPointsDispl();

        fluidPatchPointsDispl() += aitkenRelaxationFactor_*residual();
    }


    // Philip: not used on of30 as e do not have globalFaceZones
    // Make sure that displacement on all processors is equal to one
    // calculated on master processor
    // if (Pstream::parRun())
    // {
    //     if (!Pstream::master())
    //     {
    //         fluidZonePointsDispl() = vector::zero;
    //     }

    //     //- pass to all procs
    //     reduce(fluidZonePointsDispl(), sumOp<vectorField>());

        // Philip: disable on of30
//         label globalFluidZoneIndex =
//             findIndex(fluid().globalFaceZones(), fluidZoneIndex());

//         if (globalFluidZoneIndex == -1)
//         {
//             FatalErrorIn
//             (
//                 "fluidSolidInterface::updateDisplacement()"
//             )   << "global zone point map is not availabel"
//                 << abort(FatalError);
//         }

//         const labelList& map =
//             fluid().globalToLocalFaceZonePointMap()[globalFluidZoneIndex];

//         if (!Pstream::master())
//         {
//             vectorField fluidZonePointsDisplGlobal =
//                 fluidZonePointsDispl();

//             forAll(fluidZonePointsDisplGlobal, globalPointI)
//             {
//                 label localPoint = map[globalPointI];

//                 fluidZonePointsDispl()[localPoint] =
//                     fluidZonePointsDisplGlobal[globalPointI];
//             }
//         }
        //    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolidInterfaces

} // End namespace Foam

// ************************************************************************* //
