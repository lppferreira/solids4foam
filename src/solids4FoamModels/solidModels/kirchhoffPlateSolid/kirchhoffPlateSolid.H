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
    kirchhoffPlateSolid

Description
    Kirchhoff plate solver.

    Based on "rotation-free formulation" of Torlak (2006).

    Work in progress.

Author
    Philip Cardiff, UCD.  All rights reserved.

SourceFiles
    kirchhoffPlateSolid.C

\*---------------------------------------------------------------------------*/

#ifndef kirchhoffPlateSolid_H
#define kirchhoffPlateSolid_H

#include "solidModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "uniformDimensionedFields.H"
#include "areaFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class kirchhoffPlateSolid Declaration
\*---------------------------------------------------------------------------*/

class kirchhoffPlateSolid
:
    public solidModel
{
    // Private data

        //- Finite Area mesh
        faMesh aMesh_;

        //- Transverse displacement
        //  displacement prependicular to the plane of the plate
        areaScalarField w_;

        //- Transverse displacement volume field for visualisation
        volScalarField wVf_;

        //- Moment sum
        areaScalarField M_;

        //- Moment sum volume field for visualisation
        volScalarField MVf_;

        //- Net pressure field on the plate (perpendicular to the plate normal)
        areaScalarField p_;

        //- Angle of rotation
        areaVectorField theta_;

        //- Angle of rotation volume field for visualisation
        volVectorField thetaVf_;

        //- Gradient of angle of rotation
        areaTensorField gradTheta_;

        //- Plate density
        dimensionedScalar rho_;

        //- Plate Young's modulus
        dimensionedScalar E_;

        //- Plate Poisson's ratio
        dimensionedScalar nu_;

        //- Plate thickness
        const dimensionedScalar h_;

        //- Plate bending stiffness
        dimensionedScalar bendingStiffness_;

        //- Implicit stiffness; coefficient of the Laplacian term
        //  The value of this term only affects convergence and not the answer
        //  In the case of Hookean elasticity, (2*mu + lambda) is the optimum
        const volScalarField impK_;

        //- Implicit stiffness surface field
        surfaceScalarField impKf_;

        //- Reciprocal of impK
        //  We divide by impK every time tractionBoundarySnGrad is called so we
        //  will store the reciprocal for effficiency
        const volScalarField rImpK_;

    // Private Member Functions

        //- Check if the momentum equation is converged
        bool converged
        (
            const int iCorr,
            const lduSolverPerformance& solverPerfM,
            const lduSolverPerformance& solverPerfw,
            const areaScalarField& M,
            const areaScalarField& w
        );

        //- Disallow default bitwise copy construct
        kirchhoffPlateSolid(const kirchhoffPlateSolid&);

        //- Disallow default bitwise assignment
        void operator=(const kirchhoffPlateSolid&);


protected:

    // Protected member functions

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::LINEAR_GEOMETRY;
        }


public:

    //- Runtime type information
    TypeName("kirchhoffPlate");

    // Constructors

        //- Construct from components
        kirchhoffPlateSolid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );

    // Destructor

        virtual ~kirchhoffPlateSolid()
        {}


    // Member Functions

        // Access

            //- Each solidModel must indicate whether D or DD is the primary
            //  solution variable
            virtual volVectorField& solutionD()
            {
                notImplemented(type());

                // Keep compiler happy
                return D();
            }

            //- Plate bending stiffness
            const dimensionedScalar& bendingStiffness() const
            {
                return bendingStiffness_;
            }

            //- Plate Poisson's ration
            const dimensionedScalar& nu() const
            {
                return nu_;
            }


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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
