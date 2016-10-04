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

#include "oneWayFsiFluid.H"
#include "volFields.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(oneWayFsiFluid, 0);
addToRunTimeSelectionTable(fluidModel, oneWayFsiFluid, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

oneWayFsiFluid::oneWayFsiFluid(const fvMesh& mesh)
:
    fluidModel(this->typeName, mesh),
    U_
    (
        IOobject
        (
            "U",
            runTime().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    ),
    pMesh_(mesh),
    pointD_
    (
        IOobject
        (
            "pointD",
            runTime().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pMesh_,
        dimensionedVector("0", dimLength, vector::zero)
    ),
    p_
    (
        IOobject
        (
            "p",
            runTime().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    ),
    gradp_(fvc::grad(p_)),
    phi_
    (
        IOobject
        (
            "phi",
            runTime().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(U_) & mesh.Sf()
    ),
    transportProperties_
    (
        IOobject
        (
            "transportProperties",
            runTime().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    nu_(transportProperties_.lookup("nu")),
    rho_(transportProperties_.lookup("rho"))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const volVectorField& oneWayFsiFluid::U() const
{
    return U_;
}


const volScalarField& oneWayFsiFluid::p() const
{
    return p_;
}


tmp<vectorField> oneWayFsiFluid::patchViscousForce(const label patchID) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

    tvF() = rho_.value()*nu_.value()*U().boundaryField()[patchID].snGrad();

//     vectorField n = mesh().boundary()[patchID].nf();
//     tvF() -= n*(n&tvF());

    return tvF;
}


tmp<scalarField> oneWayFsiFluid::patchPressureForce(const label patchID) const
{
    tmp<scalarField> tpF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    tpF() = rho_.value()*p().boundaryField()[patchID];

    return tpF;
}


tmp<vectorField> oneWayFsiFluid::faceZoneViscousForce
(
    const label zoneID,
    const label patchID
) const
{
    vectorField pVF = patchViscousForce(patchID);

    tmp<vectorField> tvF
    (
        new vectorField(mesh().faceZones()[zoneID].size(), vector::zero)
    );
    vectorField& vF = tvF();

    const label patchStart = mesh().boundaryMesh()[patchID].start();

    forAll(pVF, i)
    {
        vF[mesh().faceZones()[zoneID].whichFace(patchStart + i)] = pVF[i];
    }

    // Parallel data exchange: collect pressure field on all processors
    reduce(vF, sumOp<vectorField>());

    return tvF;
}


tmp<scalarField> oneWayFsiFluid::faceZonePressureForce
(
    const label zoneID,
    const label patchID
) const
{
    scalarField pPF = patchPressureForce(patchID);

    tmp<scalarField> tpF
    (
        new scalarField(mesh().faceZones()[zoneID].size(), 0)
    );
    scalarField& pF = tpF();

    const label patchStart = mesh().boundaryMesh()[patchID].start();

    forAll(pPF, i)
    {
        pF[mesh().faceZones()[zoneID].whichFace(patchStart + i)] = pPF[i];
    }

    // Parallel data exchange: collect pressure field on all processors
    reduce(pF, sumOp<scalarField>());

    return tpF;
}


tmp<scalarField> oneWayFsiFluid::faceZoneMuEff
(
    const label zoneID,
    const label patchID
) const
{
    tmp<scalarField> tMuEff
    (
        new scalarField
        (
            mesh().faceZones()[zoneID].size(),
            rho_.value()*nu_.value()
        )
    );

    return tMuEff;
}

tmp<surfaceScalarField> oneWayFsiFluid::phi() const
{
    return phi_;
}

tmp<vectorField>
 oneWayFsiFluid::currentFaceZonePoints(const label zoneID) const
 {
     vectorField pointDisplacement
     (
         mesh().faceZones()[zoneID]().localPoints().size(),
         vector::zero
     );

     const vectorField& pointDI = pointD_.internalField();

     label globalZoneIndex = findIndex(globalFaceZones(), zoneID);

     if (globalZoneIndex != -1)
     {
         // global face zone
         const labelList& curPointMap =
             globalToLocalFaceZonePointMap()[globalZoneIndex];

         const labelList& zoneMeshPoints =
             mesh().faceZones()[zoneID]().meshPoints();

         vectorField zonePointsDisplGlobal
         (
             zoneMeshPoints.size(),
             vector::zero
         );

         //- Inter-proc points are shared by multiple procs
         //  pointNumProc is the number of procs which a point lies on
         scalarField pointNumProcs(zoneMeshPoints.size(), 0);

         forAll(zonePointsDisplGlobal, globalPointI)
         {
             label localPoint = curPointMap[globalPointI];

             if(zoneMeshPoints[localPoint] < mesh().nPoints())
             {
                 label procPoint = zoneMeshPoints[localPoint];

                 zonePointsDisplGlobal[globalPointI] =
                     pointDI[procPoint];

                 pointNumProcs[globalPointI] = 1;
             }
         }

         if (Pstream::parRun())
         {
             reduce(zonePointsDisplGlobal, sumOp<vectorField>());
             reduce(pointNumProcs, sumOp<scalarField>());

             //- now average the displacement between all procs
             zonePointsDisplGlobal /= pointNumProcs;
         }

         forAll(pointDisplacement, globalPointI)
         {
             label localPoint = curPointMap[globalPointI];

             pointDisplacement[localPoint] =
                 zonePointsDisplGlobal[globalPointI];
         }
     }
     else
     {
         pointDisplacement =
             vectorField
             (
                 pointDI,
                 mesh().faceZones()[zoneID]().meshPoints()
             );
     }

     tmp<vectorField> tCurrentPoints
     (
         new vectorField
         (
             mesh().faceZones()[zoneID]().localPoints()
           + pointDisplacement
         )
     );

     return tCurrentPoints;
 }

void oneWayFsiFluid::evolve()
{
    Info << "Evolving fluid model" << endl;

    // Read the latest fluid mesh
    const_cast<fvMesh&>(mesh()).readUpdate();

    // Read the latest velocity and pressure fields

    IOobject Uheader
    (
        "U",
        runTime().timeName(),
        mesh(),
        IOobject::MUST_READ
    );

    IOobject pheader
    (
        "p",
        runTime().timeName(),
        mesh(),
        IOobject::MUST_READ
    );

    if (Uheader.headerOk() && pheader.headerOk())
    {
        U_ = volVectorField(Uheader, mesh());

        p_ = volScalarField(pheader, mesh());
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels

} // End namespace Foam

// ************************************************************************* //
