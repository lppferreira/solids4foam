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

\*---------------------------------------------------------------------------*/

#include "mixedTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mixedTemperatureFvPatchScalarField::mixedTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    fieldName_("undefined")
{
    valueFraction() = 0.0;
    refValue() = 0.0;
    refGrad() = 0.0;
}


mixedTemperatureFvPatchScalarField::mixedTemperatureFvPatchScalarField
(
    const mixedTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    fieldName_(ptf.fieldName_)
{}


mixedTemperatureFvPatchScalarField::mixedTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict),
    fieldName_(dimensionedInternalField().name())
{}


mixedTemperatureFvPatchScalarField::mixedTemperatureFvPatchScalarField
(
    const mixedTemperatureFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf),
    fieldName_(tppsf.fieldName_)
{}


mixedTemperatureFvPatchScalarField::mixedTemperatureFvPatchScalarField
(
    const mixedTemperatureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF),
    fieldName_(tppsf.fieldName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void mixedTemperatureFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    // lookup: grad field from solver
    const fvPatchField<vector>& gradField =
        patch().lookupPatchField<volVectorField, vector>
        (
            "grad(" + fieldName_ + ")"
        );

    // Unit normals
    const vectorField n = patch().nf();

    // Delta vectors
    const vectorField delta = patch().delta();

    // non-Orthogonal correction vectors
    const vectorField k = (I - sqr(n)) & delta;

    Field<scalar>::operator=
    (
        valueFraction()*refValue()
      + (1.0 - valueFraction())*
        (
            this->patchInternalField()
          + (k & gradField.patchInternalField())
          + refGrad()/this->patch().deltaCoeffs()
        )
    );

    fvPatchField<scalar>::evaluate();
}


tmp<Field<scalar> > mixedTemperatureFvPatchScalarField::snGrad() const
{
    // lookup: grad field from solver
    const fvPatchField<vector>& gradField =
        patch().lookupPatchField<volVectorField, vector>
        (
            "grad(" + fieldName_ + ")"
        );

    // Unit normals
    const vectorField n = patch().nf();

    // Delta vectors
    const vectorField delta = patch().delta();

    // non-Orthogonal correction vectors
    const vectorField k = (I - sqr(n)) & delta;

    return
        valueFraction()
       *(
            refValue()
          - (
                this->patchInternalField()
              + (k & gradField.patchInternalField())
            )
        )
       *this->patch().deltaCoeffs()
      + (1.0 - valueFraction())*refGrad();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, mixedTemperatureFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
