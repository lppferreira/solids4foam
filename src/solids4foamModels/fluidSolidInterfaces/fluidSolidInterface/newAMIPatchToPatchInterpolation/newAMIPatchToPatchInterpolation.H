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
    newAMIPatchToPatchInterpolation

Description
    AMIPatchToPatchInterpolation including point-to-point interpolation
    functions.

    The point-to-point interpolation is achieved using
    PrimitivePatchInterpolation to interpolate from points-to-faces and
    then using the AMI to interpolate from faces-to-faces and then using
    a second PrimitivePatchInterpolation to interpolate from faces-to-points.

Author
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    newAMIPatchToPatchInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef newAMIPatchToPatchInterpolation_H
#define newAMIPatchToPatchInterpolation_H

#include "AMIPatchToPatchInterpolation.H"
#include "PrimitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class newAMIPatchToPatchInterpolation Declaration
\*---------------------------------------------------------------------------*/

class newAMIPatchToPatchInterpolation
:
    public AMIPatchToPatchInterpolation 
{
    // Private data
        
        //- Store reference to the source patch
        const PrimitivePatch<face, SubList, const pointField&>& srcPatch_;

        //- Store reference to the target patch
        const PrimitivePatch<face, SubList, const pointField&>& tgtPatch_;

        //- Patch interpolate for source patch
        mutable PrimitivePatchInterpolation
        <
            PrimitivePatch<face, SubList, const pointField&>
        >* srcPatchInterpPtr_;

        //- Patch interpolate for target patch
        mutable PrimitivePatchInterpolation
        <
            PrimitivePatch<face, SubList, const pointField&>
        >* tgtPatchInterpPtr_;


    // Private Member Functions

        //- Make the source patch interpolator
        void makeSrcPatchInterp() const;

        //- Make the target patch interpolator
        void makeTgtPatchInterp() const;

        //- Return the source patch interpolator
        const PrimitivePatchInterpolation
        <
            PrimitivePatch<face, SubList, const pointField&>
        >& srcPatchInterp() const;

        //- Return the target patch interpolator
        const PrimitivePatchInterpolation
        <
            PrimitivePatch<face, SubList, const pointField&>
        >& tgtPatchInterp() const;

        //- Delete demand driven data
        void clearout() const;

        //- Disallow default bitwise copy construct
        newAMIPatchToPatchInterpolation(const newAMIPatchToPatchInterpolation&);

        //- Disallow default bitwise assignment
        void operator=(const newAMIPatchToPatchInterpolation&);


public:

    //- Runtime type information
    TypeName("newAMIPatchToPatchInterpolation");


    // Constructors

        //- Construct from patches
        newAMIPatchToPatchInterpolation
        (
            const PrimitivePatch<face, SubList, const pointField&>& srcPatch,
            const PrimitivePatch<face, SubList, const pointField&>& tgtPatch,
            const faceAreaIntersect::triangulationMode& triMode,
            const bool requireMatch = true,
            const interpolationMethod& method = imFaceAreaWeight,
            const scalar lowWeightCorrection = -1,
            const bool reverseTarget = false
        );


    // Destructor

        virtual ~newAMIPatchToPatchInterpolation()
        {}


    // Member Functions

        // Evaluation

            //- Interpolate points from target to source
            template<class Type>
            tmp<Field<Type> > pointInterpolateToSource
            (
                const Field<Type>& fld
            ) const;

            //- Interpolate points from source to target
            template<class Type>
            tmp<Field<Type> > pointInterpolateToTarget
            (
                const Field<Type>& fld
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "newAMIPatchToPatchInterpolationTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
