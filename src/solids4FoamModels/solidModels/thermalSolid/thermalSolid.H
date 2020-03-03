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
    Mathematical model to solve the heat equation in solid region.
    Note: mesh motion is not supported.

Author
    Danial Khazaei, NIT.  All rights reserved.

SourceFiles
    thermalSolid.C

\*---------------------------------------------------------------------------*/

#ifndef thermalSolid_H
#define thermalSolid_H

#include "solidModel.H"
#include "thermalModel.H"
#ifdef OPENFOAMESIORFOUNDATION
    #include "SolverPerformance.H"
    #include "fvOptions.H"
#endif

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

        //- Solid temperature field [K]
        volScalarField T_;

        //- Gradient of the temperature [K/m]
        volVectorField gradT_;

        //- Solid density [kg/m3]
        volScalarField rho_;

        //- Solid specific heat [m2/s2-K]
        volScalarField C_;

        //- Solid thermal conductivity [kg-m/s3-K]
        volScalarField kappa_;

        //- Absolute temperature tolerance (in degrees)
        const scalar absTTol_;

#ifdef OPENFOAMESIORFOUNDATION
        //- Create finite volume options
        fv::options& fvOptions_;
#endif


    // Private Member Functions

        //- Calculate solid region Diffusion number
        void solidRegionDiffNo();

        //- Check if the momentum equation is converged
        bool converged
        (
            const int iCorr,
#ifdef OPENFOAMESIORFOUNDATION
            const SolverPerformance<scalar>& solverPerfT,
#else
            const lduSolverPerformance& solverPerfT,
#endif
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
                notImplemented("thermalSolid::solutionD()");
                return const_cast<volVectorField&>(volVectorField::null());
            }

            //- Return const temperature field
            virtual const volScalarField& T() const
            {
                return T_;
            }

            //- Return temperature field
            virtual volScalarField& T()
            {
                return T_;
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
            ) const
            {
                notImplemented("thermalSolid::tractionBoundarySnGrad() const");
                return tmp<vectorField>(new vectorField());
            }

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
