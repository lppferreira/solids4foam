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
    unsThermalNonLinGeomTotalLagSolid

Description
    Mathematical model for finite strains and rotations, where a total
    Lagrangian and total displacement approach is employed.

    The heat equation and momentum equation are solved in a strongly coupled
    fashion i.e. just like in solidDisplacementFoam.

    The stress is calculated by the run-time selectable mechanical law.

    "uns" in the class name refers to the "unstructured" and indicates that
    the face tangential gradient are calculated using a face-Gauss approach.

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.
    Philip Cardiff, UCD.  All rights reserved.

Modified by
    Danial Khazaei, NIT. All rights reserved.

SourceFiles
    unsThermalNonLinGeomTotalLagSolid.C

\*---------------------------------------------------------------------------*/

#ifndef unsThermalNonLinGeomTotalLagSolid_H
#define unsThermalNonLinGeomTotalLagSolid_H

#include "solidModel.H"
#include "thermalModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
              Class unsThermalNonLinGeomTotalLagSolid Declaration
\*---------------------------------------------------------------------------*/

class unsThermalNonLinGeomTotalLagSolid
:
    public solidModel
{
    // Private data

        //- Cauchy (true) stress surface field
        surfaceSymmTensorField sigmaf_;

        //- Surface displacement gradient
        surfaceTensorField gradDf_;

        //- Total deformation gradient
        volTensorField F_;

        //- Total deformation gradient surface field
        surfaceTensorField Ff_;

        //- Inverse of the total deformation gradient
        volTensorField Finv_;

        //- Inverse of the total deformation gradient surface field
        surfaceTensorField Finvf_;

        //- Jacobian of total deformation gradient
        volScalarField J_;

        //- Jacobian of total deformation gradient surface field
        surfaceScalarField Jf_;

        //- Thermal model
        thermalModel thermal_;

        //- The product of density and specific heat
        volScalarField rhoC_;

        //- Thermal conductivity
        volScalarField kappa_;

        //- Absolute temperature tolerance (in degrees)
        const scalar absTTol_;

        //- Implicit stiffness; coefficient of the Laplacian term
        //  The value of this term only affects convergence and not the answer
        //  In the case of Hookean elasticity, (2*mu + lambda) is the optimum
        const volScalarField impK_;

        //- Implicit stiffness surface field
        const surfaceScalarField impKf_;

        //- Reciprocal of impK
        //  We divide by impK every time tractionBoundarySnGrad is called so we
        //  will store the reciprocal for effficiency
        const volScalarField rImpK_;

        //- Switch to enable/disable nonlinear geometry
        const Switch nonLinear_;

        //- Debug switch for the print-out of residuals
        const Switch debug_;

        //- Viscosity/damping coefficient
        const dimensionedScalar K_;

        //- Relative tolerance
        const scalar relativeTol_;


    // Private Member Functions

        //- Calculate Diffusion number
        void DiffusionNo();

        //- Check if the momentum equation is converged
        bool converged
        (
            const int iCorr,
            const lduSolverPerformance& solverPerfD,
            const lduSolverPerformance& solverPerfT,
            const volVectorField& D,
            const volScalarField& T
        );

        //- Disallow default bitwise copy construct
        unsThermalNonLinGeomTotalLagSolid(const unsThermalNonLinGeomTotalLagSolid&);

        //- Disallow default bitwise assignment
        void operator=(const unsThermalNonLinGeomTotalLagSolid&);


protected:

    // Protected member functions

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::TOTAL_LAGRANGIAN;
        }


public:

    //- Runtime type information
    TypeName("unsThermalNonLinearGeometryTotalLagrangian");

    // Constructors

        //- Construct from components
        unsThermalNonLinGeomTotalLagSolid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~unsThermalNonLinGeomTotalLagSolid()
        {}


    // Member Functions

        // Access

            //- Each solidModel must indicate whether D or DD is the primary
            //  solution variable
            virtual volVectorField& solutionD()
            {
                // This model solves for D
                return D();
            }

            //- Return nonLinear geometry switch
            const Switch& nonLinear() const
            {
                return nonLinear_;
            }

            //- Deformed patch face area vector
            virtual tmp<vectorField> currentPatchSf
            (
                const label patchID
            ) const
            {
                return
                (
                    Jf_.boundaryField()[patchID]
                  * Finvf_.boundaryField()[patchID].T()
                  & mesh().boundary()[patchID].Sf()
                );
            }

            //- Patch heat flux (W/m2)
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
