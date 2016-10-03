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


// * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::newAMIPatchToPatchInterpolation::pointInterpolateToSource
(
    const Field<Type>& tgtFld
) const
{
    // Interpolate from target points to target faces
    const Field<Type> tgtFldAtFaces =
        tgtPatchInterp().pointToFaceInterpolate(tgtFld);

    // Interpolate from target faces to source faces
    const Field<Type> srcFldFaces = interpolateToSource(tgtFldAtFaces);

    // Interpolate from source faces to source points
    return
        srcPatchInterp().faceToPointInterpolate
        (
            srcFldFaces
        );
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::newAMIPatchToPatchInterpolation::pointInterpolateToTarget
(
    const Field<Type>& srcFld
) const
{
    // Interpolate from source points to source faces
    const Field<Type> srcFldAtFaces =
        srcPatchInterp().pointToFaceInterpolate(srcFld);

    // Interpolate from source faces to target faces
    const Field<Type> tgtFldFaces = interpolateToTarget(srcFldAtFaces);

    // Interpolate from target faces to target points
    return
        tgtPatchInterp().faceToPointInterpolate
        (
            tgtFldFaces
        );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
