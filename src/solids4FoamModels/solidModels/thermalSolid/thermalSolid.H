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
    thermalSolid

Description

    The heat transfer solver for solid models.

Author
    Danial Khazaei, NIT.  All rights reserved.

SourceFiles
    thermalSolid.C

\*---------------------------------------------------------------------------*/

#ifndef thermalSolid_H
#define thermalSolid_H

#include "solidModel.H"
#include "thermalModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class thermalSolid Declaration
\*---------------------------------------------------------------------------*/

class thermalSolid
:
    public solidModel
{
    // Private data

        //- Thermal model
        thermalModel thermal_;

        //- The product of density and specific heat
        volScalarField rhoC_;

        //- Thermal conductivity
        volScalarField k_;

        //- Absolute temperature tolerance (in degrees)
        const scalar absTTol_;

	//- Diffusion number
	scalar DiNo_;

    // Private Member Functions

        //- Check if the momentum equation is converged
        bool converged
        (
            const int iCorr,
            const lduSolverPerformance& solverPerfT,
            const volScalarField& T
        );

        //- Disallow default bitwise copy construct
        thermalSolid(const thermalSolid&);

        //- Disallow default bitwise assignment
        void operator=(const thermalSolid&);


protected:

    // Protected member functions

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::LINEAR_GEOMETRY;
        }


public:

    //- Runtime type information
    TypeName("thermalSolid");

    // Constructors

        //- Construct from components
        thermalSolid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~thermalSolid()
        {}


    // Member Functions

        // Access

            //- Each solidModel must indicate whether D or DD is the primary
            //  solution variable
            virtual volVectorField& solutionD()
            {
                // This model solves for D
                notImplemented("thermalSolid::solutionD()");
                return const_cast<volVectorField&>(volVectorField::null());
            }

            //- Return Diffusion Number
            virtual scalar& DiNo();

            //- Patch thermal flux
            virtual tmp<scalarField> patchThermalFlux
            (
                const label patchID
            ) const;

            //- Patch temperature (K)
            virtual tmp<scalarField> patchTemperature
            (
                const label patchID
            ) const;

            //- Patch KDelta
            virtual tmp<scalarField> patchKDelta
            (
                const label patchID
            ) const;

        // Edit

            //- Set temperature and flux at specified patch
            virtual void setTemperature
            (
                const label patchID,
                const scalarField& faceZoneTemperature,
                const scalarField& faceZoneKDelta
            );

            //- Evolve the solid solver and solve the mathematical model
            virtual bool evolve();

            //- Traction boundary surface normal gradient
            virtual tmp<vectorField> tractionBoundarySnGrad
            (
                const vectorField& traction,
                const scalarField& pressure,
                const fvPatch& patch
            ) const;

            //- Write fields
            virtual void writeFields(const Time& runTime);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
