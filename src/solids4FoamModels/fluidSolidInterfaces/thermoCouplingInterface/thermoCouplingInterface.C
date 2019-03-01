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

#include "thermoCouplingInterface.H"
#include "addToRunTimeSelectionTable.H"
#include "staticFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidSolidInterfaces
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(thermoCouplingInterface, 0);
addToRunTimeSelectionTable
(
    physicsModel, thermoCouplingInterface, fluidSolidInteraction
);
addToRunTimeSelectionTable
(
    fluidSolidInterface, thermoCouplingInterface, dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermoCouplingInterface::thermoCouplingInterface
(
    Time& runTime,
    const word& region
)
:
    fluidSolidInterface(typeName, runTime, region)
{
    if (fluidMesh().type() != staticFvMesh::typeName)
    {
        FatalErrorIn("thermoCouplingInterface(Time&, const word&)")
            << fluidMesh().type() <<  " is selected instead of "
                << staticFvMesh::typeName << "!\n"
                << "This coupling interface only supports "
                << "coupled temperature field without dynamic mesh.\n"
                << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool thermoCouplingInterface::evolve()
{
    initializeFields();

    updateInterpolatorAndGlobalPatches();

    scalar thermalResidualNorm = 0;

    do
    {
        outerCorr()++;

        // Transfer temperature and flux from the solid to the fluid
        updateFluidPatchTemperatureBC();

        // Solve fluid
        fluid().evolve();

        // Transfer temperature and flux from the fluid to the solid
        updateSolidPatchTemperatureBC();

        // Solve solid
        solid().evolve();

        // Calculate thermal residual
        thermalResidualNorm = updateThermalResidual();
    }
    while (thermalResidualNorm > thermalTolerance() && outerCorr() < nOuterCorr());

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolidInterfaces

} // End namespace Foam

// ************************************************************************* //
