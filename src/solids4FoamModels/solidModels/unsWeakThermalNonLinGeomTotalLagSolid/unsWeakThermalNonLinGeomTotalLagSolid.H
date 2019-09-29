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
    nonLinGeomTotalLagSolid

Description
    Mathematical model for finite strains and rotations, where a total
    Lagrangian and total displacement approach is employed.

    The stress is calculated by the run-time selectable mechanical law.

    The heat equation and momentum equation are solved in a weakly coupled
    fashion.

    "uns" in the class name refers to the "unstructured" and indicates that
    the face tangential gradient are calculated using a face-Gauss approach.

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.
    Philip Cardiff, UCD.  All rights reserved.

Modified by
    Danial Khazaei, NIT. All rights reserved.

SourceFiles
    unsWeakThermalNonLinGeomTotalLagSolid.C

\*---------------------------------------------------------------------------*/

#ifndef unsWeakThermalNonLinGeomTotalLagSolid_H
#define unsWeakThermalNonLinGeomTotalLagSolid_H

#include "unsNonLinGeomTotalLagSolid.H"
#include "thermalModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
           Class unsWeakThermalNonLinGeomTotalLagSolid Declaration
\*---------------------------------------------------------------------------*/

class unsWeakThermalNonLinGeomTotalLagSolid
:
    public unsNonLinGeomTotalLagSolid
{
    // Private data

        //- Interface patch in its deformed configuration
        mutable standAlonePatch* deformedPatchPtr_;

        //- Thermal model
        thermalModel thermal_;

        //- The product of density and specific heat
        volScalarField rhoC_;

        //- Thermal conductivity
        volScalarField kappa_;

        //- Absolute temperature tolerance (in degrees)
        const scalar absTTol_;

	//- Diffusion number
	scalar DiffusionNo_;


    // Private Member Functions

        //- Check convergence
        bool converged
        (
            const int iCorr,
            const lduSolverPerformance& solverPerfT,
            const volScalarField& T
        );

        //- Patch in deformed configuration
        const standAlonePatch& currentBoundaryPatch
        (
            const label patchID
        ) const;

        //- Patch deltaCoeffs in deformed configuration
        tmp<scalarField> currentDeltaCoeffs
        (
            const label patchID
        ) const;

        //- Disallow default bitwise copy construct
        unsWeakThermalNonLinGeomTotalLagSolid(const unsWeakThermalNonLinGeomTotalLagSolid&);

        //- Disallow default bitwise assignment
        void operator=(const unsWeakThermalNonLinGeomTotalLagSolid&);


protected:

    // Protected member functions

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::TOTAL_LAGRANGIAN;
        }


public:

    //- Runtime type information
    TypeName("unsWeakThermalNonLinearGeometryTotalLagrangian");

    // Constructors

        //- Construct from components
        unsWeakThermalNonLinGeomTotalLagSolid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~unsWeakThermalNonLinGeomTotalLagSolid();


    // Member Functions

        // Access

            //- Patch deformed face area vector
            // Note: we have two options to calculate deformed patch face area vector
            // 1. Use the cell-centred deformation gradient field 'F'.
            // 2. Calculate the deformed face area vector by interpolating displacements
            //    to the points and calculating the face area vector on the deformed patch.
            // For now we use method 1. as the result of the method 2. is not consistent,
            // perhaps due the governing equation formulation
            virtual tmp<vectorField> patchCurrentSf
            (
                const label patchID
            ) const
            {
                return
                (
                    Jf().boundaryField()[patchID]
                  * Finvf().boundaryField()[patchID].T()
                  & mesh().boundary()[patchID].Sf()
                );

                // Note: Does not produce consistent result
                // return currentBoundaryPatch(patchID).areas();
            }

            //- Return Diffusion Number
            virtual scalar& DiffusionNo();

            //- Patch heat flux
            virtual tmp<scalarField> patchHeatFlux
            (
                const label patchID
            ) const;

            //- Patch temperature (K)
            virtual tmp<scalarField> patchTemperature
            (
                const label patchID
            ) const;

            //- Patch kappa*delta
            virtual tmp<scalarField> patchKappaDelta
            (
                const label patchID
            ) const;

        // Edit

            //- Set temperature and flux at specified patch
            virtual void setTemperature
            (
                const label patchID,
                const scalarField& nbrFaceZoneTemperature,
                const scalarField& nbrFaceZoneKappaDelta
            );

            //- Evolve the solid solver and solve the mathematical model
            virtual bool evolve();

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
