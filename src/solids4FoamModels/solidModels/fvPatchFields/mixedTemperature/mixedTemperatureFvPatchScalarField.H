/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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
    Foam::mixedTemperatureFvPatchScalarField

Description
    mixed temperature boundary condition for thermal solid
    solvers with non-Orthogonal correction.

Author
    Danial Khazaei, NIT.

SourceFiles
    mixedTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef mixedTemperatureFvPatchScalarField_H
#define mixedTemperatureFvPatchScalarField_H

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

        //- Construct by mapping given mixedTemperatureFvPatchScalarField
        // onto a new patch
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

            //- Return gradient at boundary
            virtual tmp<Field<scalar> > snGrad() const;

            //- Evaluate the patch field
            virtual void evaluate
            (
#ifdef OPENFOAMESIORFOUNDATION
                const Pstream::commsTypes commsType = Pstream::commsTypes::blocking
#else
                const Pstream::commsTypes commsType = Pstream::blocking
#endif
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
