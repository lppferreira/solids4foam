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

\*---------------------------------------------------------------------------*/

#include "newAMIPatchToPatchInterpolation.H"
//#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(newAMIPatchToPatchInterpolation, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::newAMIPatchToPatchInterpolation::makeSrcPatchInterp() const
{
    if (srcPatchInterpPtr_)
    {
        FatalErrorIn
        (
            "void Foam::newAMIPatchToPatchInterpolation::makeSrcPatchInterp() "
            "const"
        )   << "pointer already set" << abort(FatalError);
    }

    srcPatchInterpPtr_ =
        new PrimitivePatchInterpolation
        <
            PrimitivePatch<face, SubList, const pointField&>
        >(srcPatch_);
}


void Foam::newAMIPatchToPatchInterpolation::makeTgtPatchInterp() const
{
    if (tgtPatchInterpPtr_)
    {
        FatalErrorIn
        (
            "void Foam::newAMIPatchToPatchInterpolation::makeTgtPatchInterp()"
            " const"
        )   << "pointer already set" << abort(FatalError);
    }

    tgtPatchInterpPtr_ =
        new PrimitivePatchInterpolation
        <
            PrimitivePatch<face, SubList, const pointField&>
        >(tgtPatch_);
}


const Foam::PrimitivePatchInterpolation
<
    Foam::PrimitivePatch<Foam::face, Foam::SubList, const Foam::pointField&>
>& Foam::newAMIPatchToPatchInterpolation::srcPatchInterp() const
{
    if (!srcPatchInterpPtr_)
    {
        makeSrcPatchInterp();
    }

    return *srcPatchInterpPtr_;
}


const Foam::PrimitivePatchInterpolation
<
    Foam::PrimitivePatch<Foam::face, Foam::SubList, const Foam::pointField&>
>& Foam::newAMIPatchToPatchInterpolation::tgtPatchInterp() const
{
    if (!tgtPatchInterpPtr_)
    {
        makeTgtPatchInterp();
    }

    return *tgtPatchInterpPtr_;
}


void Foam::newAMIPatchToPatchInterpolation::clearout() const
{
    deleteDemandDrivenData(srcPatchInterpPtr_);
    deleteDemandDrivenData(tgtPatchInterpPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::newAMIPatchToPatchInterpolation::newAMIPatchToPatchInterpolation
(
    const PrimitivePatch<face, SubList, const pointField&>& srcPatch,
    const PrimitivePatch<face, SubList, const pointField&>& tgtPatch,
    const faceAreaIntersect::triangulationMode& triMode,
    const bool requireMatch,
    const interpolationMethod& method,
    const scalar lowWeightCorrection,
    const bool reverseTarget
)
:
    AMIPatchToPatchInterpolation
    (
        srcPatch,
        tgtPatch,
        triMode,
        requireMatch,
        method,
        lowWeightCorrection,
        reverseTarget
    ),
    srcPatch_(srcPatch),
    tgtPatch_(tgtPatch),
    srcPatchInterpPtr_(NULL),
    tgtPatchInterpPtr_(NULL)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
