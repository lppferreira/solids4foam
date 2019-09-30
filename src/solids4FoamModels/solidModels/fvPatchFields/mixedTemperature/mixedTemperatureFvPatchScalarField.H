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
    mixedTemperatureFvPatchScalarField

Description
    mixed temperature boundary condition for thermal-solid
    foam solvers with non-Orthogonal correction and support
    for Total Lagrangian and Updated Lagrangian approaches.

Author
   Danial Khazaei NIT

SourceFiles
    mixedTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef mixedTemperatureFvPatchScalarField_H
#define mixedTemperatureFvPatchScalarField_H

#include "fvPatchFields.H"
#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class mixedTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class mixedTemperatureFvPatchScalarField
:
    public mixedFvPatchScalarField
{

    // Private Data

        //- Name of the temperature field
        word fieldName_;

        //- non-Orthogonal correction vector
        vectorField nonOrthCorrVector_;

        //- Patch deltaCoeffs: deformed configuration
        scalarField deltaCoeffsCurrent_;


    // Private Member Functions

        //- Is a static mesh (e.g. total Lagrangian) approach being used
        bool staticMesh() const;

        //- Is a moving mesh (e.g. updated Lagrangian) approach being used
        bool movingMesh() const;

        //- Is a nonlinear geometry (i.e. finite strain) approach being used
        bool nonLinearGeometry() const;

        //- Calculate and update the deformed fields
        void updateDeformedFields();


public:

    //- Runtime type information
    TypeName("mixedTemperature");


    // Constructors

        //- Construct from patch and internal field
        mixedTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mixedTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given
        //  mixedTemperatureFvPatchScalarField onto a new patch
        mixedTemperatureFvPatchScalarField
        (
            const mixedTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        mixedTemperatureFvPatchScalarField
        (
            const mixedTemperatureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new mixedTemperatureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        mixedTemperatureFvPatchScalarField
        (
            const mixedTemperatureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new mixedTemperatureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType = Pstream::blocking
            );

            //- Return gradient at boundary
            virtual tmp<Field<scalar> > snGrad() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<scalar> > valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<scalar> > gradientInternalCoeffs() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<scalar> > gradientBoundaryCoeffs() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
