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
    nonLinGeomUpdatedLagSolid

Description
    Mathematical model for finite strains and rotations, where an updated
    Lagrangian and incremental displacement approach is employed.

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
    unsThermalNonLinGeomUpdatedLagSolid.C

\*---------------------------------------------------------------------------*/

#ifndef unsThermalNonLinGeomUpdatedLagSolid_H
#define unsThermalNonLinGeomUpdatedLagSolid_H

#include "solidModel.H"
#include "thermalModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
             Class unsThermalNonLinGeomUpdatedLagSolid Declaration
\*---------------------------------------------------------------------------*/

class unsThermalNonLinGeomUpdatedLagSolid
:
    public solidModel
{
    // Private data

        //- Cauchy (true) stress surface field
        surfaceSymmTensorField sigmaf_;

        //- Surface displacement increment gradient
        surfaceTensorField gradDDf_;

        //- Total deformation gradient
        volTensorField F_;

        //- Total deformation gradient surface field
        surfaceTensorField Ff_;

        //- Jacobian of total deformation gradient
        volScalarField J_;

        //- Jacobian of total deformation gradient surface field
        surfaceScalarField Jf_;

        //- Relative deformation gradient
        volTensorField relF_;

        //- Relative deformation gradient surface field
        surfaceTensorField relFf_;

        //- Relative Jacobian: Jacobian of the relative deformation gradient
        volTensorField relFinv_;

        //- Relative Jacobian: surface field
        surfaceTensorField relFinvf_;

        //- Relative Jacobian: Jacobian of the relative deformation gradient
        volScalarField relJ_;

        //- Relative Jacobian: surface field
        surfaceScalarField relJf_;

        //- Density (updated configuration)
        volScalarField rho_;

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

        //- Stabilise hydrostatic pressure
        //  Apply smoothing to the relative Jacobian calculation
        const Switch stabilisePressure_;

	//- Diffusion number
	scalar DiffusionNo_;


    // Private Member Functions

        //- Check if the momentum equation is converged
        bool converged
        (
            const int iCorr,
            const lduSolverPerformance& solverPerfDD,
            const lduSolverPerformance& solverPerfT,
            const volVectorField& DD,
            const volScalarField& T
        );

        //- Disallow default bitwise copy construct
        unsThermalNonLinGeomUpdatedLagSolid(const unsThermalNonLinGeomUpdatedLagSolid&);

        //- Disallow default bitwise assignment
        void operator=(const unsThermalNonLinGeomUpdatedLagSolid&);


protected:

    // Protected member functions

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::UPDATED_LAGRANGIAN;
        }


public:

    //- Runtime type information
    TypeName("unsThermalNonLinearGeometryUpdatedLagrangian");

    // Constructors

        //- Construct from components
        unsThermalNonLinGeomUpdatedLagSolid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~unsThermalNonLinGeomUpdatedLagSolid()
        {}


    // Member Functions

        // Access

            //- Each solidModel must indicate whether D or DD is the primary
            //  solution variable
            virtual volVectorField& solutionD()
            {
                // This model solves for DD
                return DD();
            }

            //- This solid model solves for DD instead of D
            virtual bool incremental() const
            {
                return true;
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

            //- Traction boundary surface normal gradient
            virtual tmp<vectorField> tractionBoundarySnGrad
            (
                const vectorField& traction,
                const scalarField& pressure,
                const fvPatch& patch
            ) const;

            //- Update total accumulated fields
            virtual void updateTotalFields();

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
