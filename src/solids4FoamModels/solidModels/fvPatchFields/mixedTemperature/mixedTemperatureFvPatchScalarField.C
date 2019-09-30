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

\*---------------------------------------------------------------------------*/

#include "mixedTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * * //

bool mixedTemperatureFvPatchScalarField::staticMesh() const
{
    // If the deformation gradient "F" and the displacement increment DD" are
    // found then we can assume it is a moving mesh (updated Lagrangian) case
    if (db().foundObject<volVectorField>("D") && nonLinearGeometry())
    {
        return true;
    }
    else
    {
        return false;
    }
}


bool mixedTemperatureFvPatchScalarField::movingMesh() const
{
    // If the deformation gradient "F" and the displacement increment DD" are
    // found then we can assume it is a moving mesh (updated Lagrangian) case
    if (db().foundObject<volVectorField>("DD") && nonLinearGeometry())
    {
        return true;
    }
    else
    {
        return false;
    }
}


bool mixedTemperatureFvPatchScalarField::nonLinearGeometry() const
{
    // If the deformation gradient "F" is found then we will assume the case to
    // be nonlinear geometry (i.e. finite strain)
    if (db().foundObject<volTensorField>("F"))
    {
        return true;
    }
    else
    {
        return false;
    }
}


void mixedTemperatureFvPatchScalarField::updateDeformedFields()
{
    const vectorField n = patch().nf();

    const vectorField delta = patch().delta();

    // non-Orthogonal correction vector
    // Note: using undeformed unit normal and delta vectors 
    nonOrthCorrVector_ = delta - n*(n & delta);

    if (staticMesh())
    {
        // Lookup displacement (D) from solver
        const fvPatchField<vector>& DBf =
            patch().lookupPatchField<volVectorField, vector>("D");

        // Lookup inverse of the total deformation gradient from solver
        const fvsPatchField<tensor>& FinvBf =
            patch().lookupPatchField<surfaceTensorField, tensor>("Finvf");

        // Lookup Jacobian of total deformation gradient from solver
        const fvsPatchField<scalar>& JBf =
            patch().lookupPatchField<surfaceScalarField, scalar>("Jf");

        // Unit normals: deformed configuration
        const vectorField nCurrent = JBf*FinvBf.T() & n;

        // Patch delta: deformed configuration
        const vectorField deltaCurrent =
            (DBf - DBf.patchInternalField()) + delta;

        // Patch deltaCoeffs: deformed configuration
        forAll(deltaCoeffsCurrent_, faceI)
        {
            deltaCoeffsCurrent_[faceI] =
                scalar(1) / max
                (
                    nCurrent[faceI] & deltaCurrent[faceI],
                    0.05*mag(deltaCurrent[faceI])
                );
        }
    }
    else if (movingMesh())
    {
        // Lookup displacement increment (DD) from solver
        const fvPatchField<vector>& DDBf =
            patch().lookupPatchField<volVectorField, vector>("DD");

        // Lookup inverse of the total deformation gradient from solver
        const fvsPatchField<tensor>& relFinvBf =
            patch().lookupPatchField<surfaceTensorField, tensor>("relFinvf");

        // Lookup Jacobian of total deformation gradient from solver
        const fvsPatchField<scalar>& relJBf =
            patch().lookupPatchField<surfaceScalarField, scalar>("relJf");

        // Unit normals: deformed configuration
        const vectorField nCurrent = relJBf*relFinvBf.T() & n;

        // Patch delta: deformed configuration
        const vectorField deltaCurrent =
            (DDBf - DDBf.patchInternalField()) + delta;

        // Patch deltaCoeffs: deformed configuration
        forAll(deltaCoeffsCurrent_, faceI)
        {
            deltaCoeffsCurrent_[faceI] =
                scalar(1) / max
                (
                    nCurrent[faceI] & deltaCurrent[faceI],
                    0.05*mag(deltaCurrent[faceI])
                );
        }
    }
    else
    {
        deltaCoeffsCurrent_ = this->patch().deltaCoeffs();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mixedTemperatureFvPatchScalarField::mixedTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    fieldName_("undefined"),
    nonOrthCorrVector_(p.size(), vector::zero),
    deltaCoeffsCurrent_(p.size(), 0)
{}


mixedTemperatureFvPatchScalarField::mixedTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict),
    fieldName_(dimensionedInternalField().name()),
    nonOrthCorrVector_(p.size(), vector::zero),
    deltaCoeffsCurrent_(p.size(), 0)
{
    // Call evaluate only if the value is not found. Used to avoid evaluating
    // when we have incomplete meshes during Parallel Load Balancing. When
    // shipping the field over to another processor, we first call write, making
    // sure that the value is written and read it on the other side (see
    // write member function). If this proves to be problematic, we can always
    // initialize with patch internal field for the start-up. VV, 12/Apr/2019.
    if (!dict.found("value"))
    {
        evaluate();
    }
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
    fieldName_(ptf.fieldName_),
    nonOrthCorrVector_(ptf.nonOrthCorrVector_, mapper),
    deltaCoeffsCurrent_(ptf.deltaCoeffsCurrent_, mapper)
{}


mixedTemperatureFvPatchScalarField::mixedTemperatureFvPatchScalarField
(
    const mixedTemperatureFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    fieldName_(ptf.fieldName_),
    nonOrthCorrVector_(ptf.nonOrthCorrVector_),
    deltaCoeffsCurrent_(ptf.deltaCoeffsCurrent_)
{}


mixedTemperatureFvPatchScalarField::mixedTemperatureFvPatchScalarField
(
    const mixedTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    fieldName_(ptf.fieldName_),
    nonOrthCorrVector_(ptf.nonOrthCorrVector_),
    deltaCoeffsCurrent_(ptf.deltaCoeffsCurrent_)
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

    // Update deformed fields
    updateDeformedFields();

    // Lookup grad from solver
    const fvPatchField<vector>& gradField =
        patch().lookupPatchField<volVectorField, vector>
        (
            "grad(" + fieldName_ + ")"
        );

    Field<scalar>::operator=
    (
        valueFraction()*refValue()
      + (1.0 - valueFraction())*
        (
            this->patchInternalField()
          + (nonOrthCorrVector_ & gradField.patchInternalField())
          + refGrad()/deltaCoeffsCurrent_
        )
    );

    fvPatchField<scalar>::evaluate();
}


tmp<Field<scalar> > mixedTemperatureFvPatchScalarField::snGrad() const
{
    // Lookup grad from solver
    const fvPatchField<vector>& gradField =
        patch().lookupPatchField<volVectorField, vector>
        (
            "grad(" + fieldName_ + ")"
        );

    return
        valueFraction()
       *(
            refValue()
          - (
                this->patchInternalField()
              + (nonOrthCorrVector_ & gradField.patchInternalField())
            )
        )
       *deltaCoeffsCurrent_
      + (1.0 - valueFraction())*refGrad();
}


tmp<Field<scalar> > mixedTemperatureFvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return
         valueFraction()*refValue()
       + (1.0 - valueFraction())*refGrad()/deltaCoeffsCurrent_;
}


tmp<Field<scalar> > mixedTemperatureFvPatchScalarField::gradientInternalCoeffs() const
{
    return -pTraits<scalar>::one*valueFraction()*deltaCoeffsCurrent_;
}


tmp<Field<scalar> > mixedTemperatureFvPatchScalarField::gradientBoundaryCoeffs() const
{
    return
        valueFraction()*deltaCoeffsCurrent_*refValue()
      + (1.0 - valueFraction())*refGrad();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    mixedTemperatureFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
