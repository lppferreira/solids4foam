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

Class
    buoyantBoussinesqPimpleFluid

Description
    Transient solver for buoyant, turbulent flow of incompressible fluids.

    Implements the PIMPLE (transient PISO) solution methodology.

    Uses the Boussinesq approximation:
    \f[
        rho_{eff} = 1 - beta(T - T_{ref})
    \f]

    where:
        \f$ rho_{eff} \f$ = the effective (driving) density
        beta = thermal expansion coefficient [1/K]
        T = temperature [K]
        \f$ T_{ref} \f$ = reference temperature [K]

    Valid when:
    \f[
        rho_{eff} << 1
    \f]

Author
    Philip Cardiff, UCD. All rights reserved.
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.

SourceFiles
    buoyantBoussinesqPimpleFluid.C

\*---------------------------------------------------------------------------*/

#ifndef buoyantBoussinesqPimpleFluid_H
#define buoyantBoussinesqPimpleFluid_H

#include "fluidModel.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"
#include "pimpleControl.H"
#include "fvm.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

/*---------------------------------------------------------------------------*\
                        Class buoyantBoussinesqPimpleFluid Declaration
\*---------------------------------------------------------------------------*/

class buoyantBoussinesqPimpleFluid
:
    public fluidModel
{
    // Private data

        //- Temperature
        volScalarField T_;

        //- Transport model
        singlePhaseTransportModel laminarTransport_;

        //- Thermal expansion coefficient [1/K]
        dimensionedScalar beta_;

        // Reference temperature [K]
        dimensionedScalar TRef_;

        // Laminar Prandtl number
        dimensionedScalar Pr_;

        // Turbulent Prandtl number
        dimensionedScalar Prt_;

        //- Density of phase 1
        autoPtr<incompressible::RASModel> turbulence_;

        //- Density
        dimensionedScalar rho_;

        //- Gravitional acceleration
        const uniformDimensionedVectorField g_;

        //- Beta*(g*h) field
        surfaceScalarField betaghf_;

        //- Kinematic density for buoyancy force
        volScalarField rhok_;


    // Private Member Functions

        //- Solve velocity equation
        tmp<fvVectorMatrix> solveUEqn();

        //- Solve temperature equation
        void solveTEqn();

        //- Solve pressure equation
        void solvePEqn
        (
            pimpleControl& pimple,
            tmp<fvVectorMatrix>& UEqn,
            const label pdRefCell,
            const scalar pdRefValue
        );

        //- Disallow default bitwise copy construct
        buoyantBoussinesqPimpleFluid(const buoyantBoussinesqPimpleFluid&);

        //- Disallow default bitwise assignment
        void operator=(const buoyantBoussinesqPimpleFluid&);


public:

    //- Runtime type information
    TypeName("buoyantBoussinesqPimpleFluid");

    // Constructors

        //- Construct from components
        buoyantBoussinesqPimpleFluid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~buoyantBoussinesqPimpleFluid()
        {}


    // Member Functions

        // Access

            //- Patch viscous force (N/m2)
            virtual tmp<vectorField> patchViscousForce
            (
                const label patchID
            ) const;

            //- Patch pressure force (N/m2)
            virtual tmp<scalarField> patchPressureForce
            (
                const label patchID
            ) const;

            //- Face zone effective dynamic viscosity
            virtual tmp<scalarField> faceZoneMuEff
            (
                const label zoneID,
                const label patchID
            ) const;

        // Edit

            //- Evolve the fluid model
            virtual bool evolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
