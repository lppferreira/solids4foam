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

\*----------------------------------------------------------------------------*/

#include "solidStressDisplacements.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidStressDisplacements, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        solidStressDisplacements,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::solidStressDisplacements::writeData()
{
    // Lookup the solid mesh
    const fvMesh* meshPtr = NULL;
    if (time_.foundObject<fvMesh>("solid"))
    {
        meshPtr = &(time_.lookupObject<fvMesh>("solid"));
    }
    else
    {
        meshPtr = &(time_.lookupObject<fvMesh>("region0"));
    }
    const fvMesh& mesh = *meshPtr;

    // Lookup the displacement field
    const vectorField& D =
        mesh.lookupObject<volVectorField>
        (
            "D"
        ).boundaryField()[historyPatchID_];

    // Calculate the average stress on the patch
    symmTensor averageStress = symmTensor::zero;

    // Check if it is a linear or nonlinear geometry case
    if (mesh.foundObject<volSymmTensorField>("sigmaCauchy"))
    {
        // Lookup the Cauchy stress
        const symmTensorField& sigma =
            mesh.lookupObject<volSymmTensorField>
            (
                "sigmaCauchy"
            ).boundaryField()[historyPatchID_];

        averageStress = gAverage(sigma);
    }
    else
    {
        // Linear geometry

        // Lookup the stress
        const symmTensorField& sigma =
            mesh.lookupObject<volSymmTensorField>
            (
                "sigma"
            ).boundaryField()[historyPatchID_];

        averageStress = gAverage(sigma);
    }

    // Arithmetic average disp and force on patch
    vector avDisp = average(D);

    if (Pstream::master())
    {
        historyFilePtr_()
            << time_.time().value() << " "
                << avDisp.x() << " "
                << avDisp.y() << " "
                << avDisp.z() << " "
                << averageStress.xx() << " "
                << averageStress.xy() << " "
                << averageStress.xz() << " "
                << averageStress.yy() << " "
                << averageStress.yz() << " "
                << averageStress.zz();

        historyFilePtr_() << endl;
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidStressDisplacements::solidStressDisplacements
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    historyPatchID_(-1),
    historyFilePtr_(NULL)
{
    Info << "Creating " << this->name() << " function object." << endl;

    word historyPatchName("notSpecified");
    if (dict.found("historyPatch"))
    {
        dict.lookup("historyPatch") >> historyPatchName;
    }
    else
    {
        FatalErrorIn(this->name() + " function object constructor")
            << "historyPatch not specified."
            << abort(FatalError);
    }

    // Lookup the solid mesh
    const fvMesh* meshPtr = NULL;
    if (time_.foundObject<fvMesh>("solid"))
    {
        meshPtr = &(time_.lookupObject<fvMesh>("solid"));
    }
    else
    {
        meshPtr = &(time_.lookupObject<fvMesh>("region0"));
    }
    const fvMesh& mesh = *meshPtr;

    historyPatchID_ = mesh.boundaryMesh().findPatchID(historyPatchName);

    if (historyPatchID_ == -1)
    {
        FatalErrorIn(this->name() + " function object constructor")
            << "history patch " << historyPatchName << " not found."
            << abort(FatalError);
    }

    // Create history file if not already created
    if (historyFilePtr_.empty())
    {
        // File update
        if (Pstream::master())
        {
            fileName historyDir;

            word startTimeName =
                time_.timeName(mesh.time().startTime().value());


            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                historyDir = time_.path()/".."/"history"/startTimeName;
            }
            else
            {
                historyDir = time_.path()/"history"/startTimeName;
            }

            // Create directory if does not exist.
            mkDir(historyDir);

            // Open new file at start up
            historyFilePtr_.reset
            (
                new OFstream
                (
                    historyDir/"solidStressDisplacements"
                  + historyPatchName + ".dat"
                )
            );

            // Add headers to output data
            if (historyFilePtr_.valid())
            {
                historyFilePtr_()
                    << "# Time" << " "
                    << "dispX" << " " << "dispY" << " "
                    << "dispZ" << " "
                    << "stressXX" << " " << "stressXY" << " "
                    << "stressXZ" << " "
                    << "stressYY" << " " << "stressYZ" << " "
                    << "stressZZ";

                historyFilePtr_() << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solidStressDisplacements::start()
{
    return writeData();
}


bool Foam::solidStressDisplacements::execute()
{
    return writeData();
}


bool Foam::solidStressDisplacements::read(const dictionary& dict)
{
    return true;
}

// ************************************************************************* //
