/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "timeVaryingVelocityPointPatchVectorField.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

timeVaryingVelocityPointPatchVectorField::
timeVaryingVelocityPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    timeSeries_(),
    curTimeIndex_(-1)
{}


timeVaryingVelocityPointPatchVectorField::
timeVaryingVelocityPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    timeSeries_(dict),
    curTimeIndex_(-1)
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }
}


timeVaryingVelocityPointPatchVectorField::
timeVaryingVelocityPointPatchVectorField
(
    const timeVaryingVelocityPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    timeSeries_(ptf.timeSeries_),
    curTimeIndex_(ptf.curTimeIndex_)
{}


timeVaryingVelocityPointPatchVectorField::
timeVaryingVelocityPointPatchVectorField
(
    const timeVaryingVelocityPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    timeSeries_(ptf.timeSeries_),
    curTimeIndex_(ptf.curTimeIndex_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void timeVaryingVelocityPointPatchVectorField::autoMap
(
    const pointPatchFieldMapper& m
)
{
    fixedValuePointPatchField<vector>::autoMap(m);
}


void timeVaryingVelocityPointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    const timeVaryingVelocityPointPatchVectorField& oVptf =
        refCast<const timeVaryingVelocityPointPatchVectorField>(ptf);

    fixedValuePointPatchField<vector>::rmap(oVptf, addr);
}


void timeVaryingVelocityPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    vector velocity = vector::zero;

    if (curTimeIndex_ != db().time().timeIndex())
    {
        curTimeIndex_ = db().time().timeIndex();

        // Lookup velocity from interpolation table
        velocity = timeSeries_(db().time().timeOutputValue());
    }

    Field<vector>::operator=(velocity);

    fixedValuePointPatchField<vector>::updateCoeffs();
}


void timeVaryingVelocityPointPatchVectorField::write(Ostream& os) const
{
    pointPatchField<vector>::write(os);
    timeSeries_.write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    timeVaryingVelocityPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
