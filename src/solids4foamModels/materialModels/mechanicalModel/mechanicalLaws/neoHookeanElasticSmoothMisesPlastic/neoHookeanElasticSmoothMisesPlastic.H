/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    neoHookeanElasticSmoothMisesPlastic

Description
    Hyper-elasto-plastic mechanical law with neo-Hookean elasticity and Mises/J2
    plasticity.

    A smooth stress return is implemented as developed by M. Rubin et al.

    Rubin and Attia, CALCULATION OF HYPERELASTIC RESPONSE OF FINITELY
    DEFORMED ELASTIC-VISCOPLASTIC MATERIALS, INTERNATIONAL JOURNAL FOR
    NUMERICAL METHODS IN ENGINEERING, VOL. 39,309-320(1996)

    M. Hollenstein M. Jabareen M. B. Rubin, Modeling a smooth elastic-
    inelastic transition with a strongly objective numerical integrator
    needing no iteration, Comput Mech (2013) 52:649–667
    DOI 10.1007/s00466-013-0838-7

    The actual stress calculation is performed by the Rubin et al. fortran
    subroutine, which is compiled and linked as an extern function.

SourceFiles
    neoHookeanElasticSmoothMisesPlastic.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef neoHookeanElasticSmoothMisesPlastic_H
#define neoHookeanElasticSmoothMisesPlastic_H

#include "mechanicalLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearElastic Declaration
\*---------------------------------------------------------------------------*/

class neoHookeanElasticSmoothMisesPlastic
:
    public mechanicalLaw
{
    // Private data

        //- Initial density
        const dimensionedScalar rho_;

        //- Young's modulus
        const dimensionedScalar E_;

        //- Poisson's ratio
        const dimensionedScalar nu_;

        //- Shear modulus
        const dimensionedScalar mu_;

        //- Bulk modulus
        const dimensionedScalar K_;

        //- Plasticity parameter
        const scalar a0_;

        //- Plasticity parameter
        const scalar a1_;

        //- Plasticity parameter
        const scalar b0_;

        //- Plasticity parameter
        const scalar b1_;

        //- Plasticity parameter - rate at which kappas is reached
        const scalar m_;

        //- Plasticity parameter - measure of asymptotic yield stress
        const scalar kappas_;

        // Relative deformation gradient volField
        volTensorField* relFPtr_;

        // Relative deformation gradient surfaceField
        surfaceTensorField* relFfPtr_;

        // Jacobian of the toal deformation gradient volField
        volScalarField* JPtr_;

        // Jacobian of the toal deformation gradient surfaceField
        surfaceScalarField* JfPtr_;

        //- Elastic left Cauchy-Green trial strain tensor with volumetric term
        //  removed. Trial means that plasticity is neglected.
        volSymmTensorField bEbarTrial_;

        //- Elastic left Cauchy-Green trial strain tensor with volumetric term
        //  removed surface field
        surfaceSymmTensorField bEbarTrialf_;

        //- Elastic left Cauchy-Green strain tensor with volumetric term removed
        volSymmTensorField bEbar_;

        //- Elastic left Cauchy-Green strain tensor with volumetric term removed
        //  surface field
        surfaceSymmTensorField bEbarf_;

        //- Current hardening parameter
        volScalarField kappa_;

        //- Current hardening parameter surface field
        surfaceScalarField kappaf_;

        //- Active yielding flag
        //     1.0 for active yielding
        //     0.0 otherwise
        volScalarField activeYield_;

        //- Maximum allowed error in the plastic strain integration
        //const scalar maxDeltaErr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        neoHookeanElasticSmoothMisesPlastic
        (
            const neoHookeanElasticSmoothMisesPlastic&
        );

        //- Disallow default bitwise assignment
        void operator=(const neoHookeanElasticSmoothMisesPlastic&);

        //- Make the relF field
        void makeRelF();

        //- Return a reference to the relF field
        volTensorField& relF();

        //- Make the relFf field
        void makeRelFf();

        //- Return a reference to the relFf field
        surfaceTensorField& relFf();

        //- Make the J field
        void makeJ();

        //- Return a reference to the J field
        volScalarField& J();

        //- Make the Jf field
        void makeJf();

        //- Return a reference to the Jf field
        surfaceScalarField& Jf();

        //- Calculate current stress
        void calculateStress
        (
            const tensor& relF,         // Relative deformation gradient
            const symmTensor& bEbarOld, // Old bEbar
            const scalar JOld,          // Old total Jacobian
            const scalar kappaOld,      // Old hardening parameter
            symmTensor& bEbar,          // New bEbar
            scalar& J,                  // New total Jacobian
            scalar& kappa,              // New hardening parameter
            symmTensor& sigma           // New Cauchy stress
        ) const;

        //- Test the calculate stress function
        //  Stress and specified relF written to file
        void testCalculateStress() const;


public:

    //- Runtime type information
    TypeName("neoHookeanElasticSmoothMisesPlastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        neoHookeanElasticSmoothMisesPlastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~neoHookeanElasticSmoothMisesPlastic();


    // Member Functions

        //- Return density
        virtual tmp<volScalarField> rho() const;

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

        //- Update the stress
        virtual void correct(volSymmTensorField& sigma);

        //- Update the stress surface field
        virtual void correct(surfaceSymmTensorField& sigma);

        //- Return material residual i.e. a measured of how convergence of
        //  the material model
        virtual scalar residual();

        //- Update the yield stress: called at end of time-step
        virtual void updateTotalFields();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
