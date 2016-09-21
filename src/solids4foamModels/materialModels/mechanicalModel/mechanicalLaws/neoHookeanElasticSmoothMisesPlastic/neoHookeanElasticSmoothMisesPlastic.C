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

#include "neoHookeanElasticSmoothMisesPlastic.H"
#include "addToRunTimeSelectionTable.H"
#include "transformGeometricField.H"
#include "fvc.H"


// * * * * * * * * * * * * * * Extern Functions  * * * * * * * * * * * * * * //

namespace Foam
{
    // Declare fortran function prototypes
    extern "C"
    {
        void smoothmultiphase_
        (
            const double* deltaT,
            const double* K,
            const double* mu,
            const double* a0,
            const double* b0,
            const double* a1,
            const double* b1,
            const double* m,
            const double* kappas,
            const double* JOld,
            const double bEbarOld[3][3],
            const double* kappaOld,
            const double relF[3][3],
            double* J,
            double bEbar[3][3],
            double* kappa,
            double sigmaGPa[3][3]
        );
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(neoHookeanElasticSmoothMisesPlastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, neoHookeanElasticSmoothMisesPlastic, dictionary
    );
} // End of namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::neoHookeanElasticSmoothMisesPlastic::makeRelF()
{
    if (relFPtr_)
    {
        FatalErrorIn
        (
            "void Foam::neoHookeanElasticSmoothMisesPlastic::makeRelF()"
        )   << "pointer already set" << abort(FatalError);
    }

    relFPtr_ =
        new volTensorField
        (
            IOobject
            (
                "relF",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedTensor("I", dimless, I)
        );
}


Foam::volTensorField& Foam::neoHookeanElasticSmoothMisesPlastic::relF()
{
    if (!relFPtr_)
    {
        makeRelF();
    }

    return *relFPtr_;
}


void Foam::neoHookeanElasticSmoothMisesPlastic::makeRelFf()
{
    if (relFfPtr_)
    {
        FatalErrorIn
        (
            "void Foam::neoHookeanElasticSmoothMisesPlastic::makeRelFf()"
        )   << "pointer already set" << abort(FatalError);
    }

    relFfPtr_ =
        new surfaceTensorField
        (
            IOobject
            (
                "relFf",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedTensor("I", dimless, I)
        );
}


Foam::surfaceTensorField& Foam::neoHookeanElasticSmoothMisesPlastic::relFf()
{
    if (!relFfPtr_)
    {
        makeRelFf();
    }

    return *relFfPtr_;
}


void Foam::neoHookeanElasticSmoothMisesPlastic::makeJ()
{
    if (JPtr_)
    {
        FatalErrorIn("void Foam::neoHookeanElasticSmoothMisesPlastic::makeJ()")
            << "pointer already set" << abort(FatalError);
    }

    JPtr_ =
        new volScalarField
        (
            IOobject
            (
                "J",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("one", dimless, 1.0)
        );
}


Foam::volScalarField& Foam::neoHookeanElasticSmoothMisesPlastic::J()
{
    if (!JPtr_)
    {
        makeJ();
    }

    return *JPtr_;
}


void Foam::neoHookeanElasticSmoothMisesPlastic::makeJf()
{
    if (JfPtr_)
    {
        FatalErrorIn("void Foam::neoHookeanElasticSmoothMisesPlastic::makeJf()")
            << "pointer already set" << abort(FatalError);
    }

    JfPtr_ =
        new surfaceScalarField
        (
            IOobject
            (
                "Jf",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("one", dimless, 1.0)
        );
}


Foam::surfaceScalarField& Foam::neoHookeanElasticSmoothMisesPlastic::Jf()
{
    if (!JfPtr_)
    {
        makeJf();
    }

    return *JfPtr_;
}


void Foam::neoHookeanElasticSmoothMisesPlastic::calculateStress
(
    const tensor& relF,         // Relative deformation gradient
    const symmTensor& bEbarOld, // Old bEbar
    const scalar JOld,          // Old total Jacobian
    const scalar kappaOld,      // Old hardening parameter
    symmTensor& bEbar,          // New bEbar
    scalar& J,                  // New total Jacobian
    scalar& kappa,              // New hardening parameter
    symmTensor& sigma           // New Cauchy stress
) const
{
    // Convert relF into an 2-D array to be passed to the Fortran
    // subroutine
    // Note: arrays are transposed when passed to Fortran
    const double fortranRelF[3][3] =
    {
        relF.xx(), relF.yx(), relF.zx(),
        relF.xy(), relF.yy(), relF.zy(),
        relF.xz(), relF.yz(), relF.zz()
    };
    const double fortranBEbarOld[3][3] =
    {
        bEbarOld.xx(), bEbarOld.xy(), bEbarOld.xz(),
        bEbarOld.xy(), bEbarOld.yy(), bEbarOld.yz(),
        bEbarOld.xz(), bEbarOld.yz(), bEbarOld.zz()
    };

    // Stress and bEbar tensors to be calculated
    double fortranSigmaGPa[3][3] = {{0}};
    double fortranBEbar[3][3] = {{0}};

    // We will pass a pointer to deltaT so we create it here before passing it,
    // rather than passing a copy which could cause memory errors
    const scalar deltaT = mesh().time().deltaTValue();

    // Calculate stress using Fortran stress return
    smoothmultiphase_
    (
        &deltaT,          // time-step size
        &K_.value(),      // bulk modulus
        &mu_.value(),     // shear modulus
        &a0_,             // plasticity parameter
        &b0_,             // plasticity parameter
        &a1_,             // plasticity parameter
        &b1_,             // plasticity parameter
        &m_,              // plasticity parameter
        &kappas_,         // asymptotic plasticity parameter
        &JOld,            // old total Jacobian
        fortranBEbarOld,  // old bEbar
        &kappaOld,        // old hardening parameter
        fortranRelF,      // relative deformation gradient
        &J,               // new total Jacobian
        fortranBEbar,     // new bEbar
        &kappa,           // new hardening parameter
        fortranSigmaGPa   // Cauchy stress to be calculated (in GPa)
    );

    // Copy sigma 2-D array into symmTensor
    // Remember: arrays are transposed when received from Fortran
    sigma =
        1e9
       *symmTensor
        (
            fortranSigmaGPa[0][0], fortranSigmaGPa[1][0], fortranSigmaGPa[2][0],
                                   fortranSigmaGPa[1][1], fortranSigmaGPa[2][1],
                                                          fortranSigmaGPa[2][2]
        );

    // Copy bEbar 2-D array into symmTensor
    bEbar =
        symmTensor
        (
            fortranBEbar[0][0], fortranBEbar[1][0], fortranBEbar[2][0],
                                fortranBEbar[1][1], fortranBEbar[2][1],
                                                    fortranBEbar[2][2]
        );
}


void Foam::neoHookeanElasticSmoothMisesPlastic::testCalculateStress() const
{
    Info<< "Testing stress point calculation" << endl;

    OFstream outFile("stressVsEngStrain.txt");

    tensor relF = I;
    symmTensor bEbarOld = I;
    symmTensor bEbar = bEbarOld;
    scalar JOld = 1.0;
    scalar J = JOld;
    scalar kappaOld = 0.01;
    scalar kappa = 0.01;
    symmTensor sigma = symmTensor::zero;

    scalar engStrain = 0.0;

    for (int i = 0; i < 500; i++)
    {
        //Info<< "i is " << i << endl;

        if (i < 100)
        {
            relF.yy() = 1.001;
            engStrain += 0.001;
        }
        else if (i < 300)
        {
            relF.yy() = 0.999;
            engStrain -= 0.001;
        }
        else if (i < 500)
        {
            relF.yy() = 1.001;
            engStrain += 0.001;
        }
        else
        {
            break;
        }

        // Info<< "relF " << relF << nl
        //     << "bEbarOld " << bEbarOld << nl
        //     << "JOld " << JOld << nl
        //     << "bEbar " << bEbar << nl
        //     << "J " << J << nl
        // Info<< "kappa " << kappa << endl;
        // Info<< "kappaOld " << kappaOld << endl;
        //kappaOld = 0.01;
        //     << "sigma " << sigma << nl
        //     << endl;

        calculateStress
        (
            relF,
            bEbarOld,
            JOld,
            kappaOld,
            bEbar,
            J,
            kappa,
            sigma
        );

        // Update old values
        bEbarOld = bEbar;
        JOld = J;
        kappaOld = kappa;

        outFile
            << 1e-9*sigma.yy() << " " << engStrain << endl;
        Info<< kappa << " " << kappaOld << " "
            << 1e-9*sigma.yy() << " " << engStrain << endl;
    }

    FatalErrorIn
    (
        "void Foam::neoHookeanElasticSmoothMisesPlastic::correct\n"
        "(\n"
        "volSymmTensorField&\n"
        ")"
    )   << "TESTING stress point calculation: writing stressVsEngStrain.txt"
        << abort(FatalError);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::neoHookeanElasticSmoothMisesPlastic::neoHookeanElasticSmoothMisesPlastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mechanicalLaw(name, mesh, dict),
    rho_(dict.lookup("rho")),
    E_(dict.lookup("E")),
    nu_(dict.lookup("nu")),
    mu_(E_/(2.0*(1.0 + nu_))),
    K_
    (
        planeStress()
      ? (nu_*E_/((1.0 + nu_)*(1.0 - nu_))) + (2.0/3.0)*mu_
      : (nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))) + (2.0/3.0)*mu_
    ),
    a0_(readScalar(dict.lookup("a0"))),
    a1_(readScalar(dict.lookup("a1"))),
    b0_(readScalar(dict.lookup("b0"))),
    b1_(readScalar(dict.lookup("b1"))),
    m_(readScalar(dict.lookup("m"))),
    kappas_(readScalar(dict.lookup("kappas"))),
    relFPtr_(NULL),
    relFfPtr_(NULL),
    JPtr_(NULL),
    JfPtr_(NULL),
    bEbarTrial_
    (
        IOobject
        (
            "bEbarTrial",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("I", dimless, I)
    ),
    bEbarTrialf_
    (
        IOobject
        (
            "bEbarTrialf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("I", dimless, I)
    ),
    bEbar_
    (
        IOobject
        (
            "bEbar",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("I", dimless, I)
    ),
    bEbarf_
    (
        IOobject
        (
            "bEbarf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("I", dimless, I)
    ),
    kappa_
    (
        IOobject
        (
            "kappa",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, readScalar(dict.lookup("kappa")))
    ),
    kappaf_
    (
        IOobject
        (
            "kappaf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, readScalar(dict.lookup("kappa")))
    ),
    activeYield_
    (
        IOobject
        (
            "activeYield",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0)
    )
    // maxDeltaErr_
    // (
    //    mesh.time().controlDict().lookupOrDefault<scalar>("maxDeltaErr", 0.01)
    // )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::neoHookeanElasticSmoothMisesPlastic::
~neoHookeanElasticSmoothMisesPlastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::neoHookeanElasticSmoothMisesPlastic::rho() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "rhoLaw",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            rho_,
            calculatedFvPatchScalarField::typeName
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::neoHookeanElasticSmoothMisesPlastic::impK() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "impK",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            (4.0/3.0)*mu_ + K_, // == 2*mu + lambda
            zeroGradientFvPatchScalarField::typeName
        )
    );
}


void Foam::neoHookeanElasticSmoothMisesPlastic::correct
(
    volSymmTensorField& sigma
)
{
    //testCalculateStress();

    if (mesh().foundObject<volTensorField>("grad(DD)"))
    {
        // Lookup gradient of displacement increment
        const volTensorField& gradDD =
            mesh().lookupObject<volTensorField>("grad(DD)");

        // Update the relative deformation gradient
        relF() = I + gradDD.T();
    }
    else
    {
        // Lookup gradient of displacement
        const volTensorField& gradD =
            mesh().lookupObject<volTensorField>("grad(D)");

        // Calculate the gradient of the displacement increment
        const volTensorField gradDD = gradD - gradD.oldTime();

        // Update the relative deformation gradient
        relF() = I + gradDD.T();
    }

    // Store field to calculate residual
    kappa_.storePrevIter();

    // Take references to the internal fields for efficiency
    const tensorField& relFI = relF().internalField();
    const symmTensorField& bEbarOldI = bEbar_.oldTime().internalField();
    const scalarField& JOldI = J().oldTime().internalField();
    const scalarField& kappaOldI = kappa_.oldTime().internalField();
    symmTensorField& bEbarI = bEbar_.internalField();
    scalarField& JI = J().internalField();
    scalarField& kappaI = kappa_.internalField();
    symmTensorField& sigmaI = sigma.internalField();

    // Calculate stress for each cell
    forAll(sigmaI, cellI)
    {
        calculateStress
        (
            relFI[cellI],
            bEbarOldI[cellI],
            JOldI[cellI],
            kappaOldI[cellI],
            bEbarI[cellI],
            JI[cellI],
            kappaI[cellI],
            sigmaI[cellI]
        );
    }

    // Calculate stress for each boundary face
    forAll(sigma.boundaryField(), patchI)
    {
        if
        (
            !sigma.boundaryField()[patchI].coupled()
         && sigma.boundaryField()[patchI].type() != "empty"
        )
        {
            // Take references to the boundary patch fields for efficiency
            const tensorField& relFP = relF().boundaryField()[patchI];
            const symmTensorField& bEbarOldP =
                bEbar_.oldTime().boundaryField()[patchI];
            const scalarField& JOldP = J().oldTime().boundaryField()[patchI];
            const scalarField& kappaOldP =
                kappa_.oldTime().boundaryField()[patchI];
            symmTensorField& bEbarP = bEbar_.boundaryField()[patchI];
            scalarField& JP = J().boundaryField()[patchI];
            scalarField& kappaP = kappa_.boundaryField()[patchI];
            symmTensorField& sigmaP = sigma.boundaryField()[patchI];

            forAll(sigmaP, faceI)
            {
                calculateStress
                (
                    relFP[faceI],
                    bEbarOldP[faceI],
                    JOldP[faceI],
                    kappaOldP[faceI],
                    bEbarP[faceI],
                    JP[faceI],
                    kappaP[faceI],
                    sigmaP[faceI]
                );
            }
        }
    }

    // Correct coupled boundaries
    bEbar_.correctBoundaryConditions();
    J().correctBoundaryConditions();
    kappa_.correctBoundaryConditions();
    sigma.correctBoundaryConditions();
}


void Foam::neoHookeanElasticSmoothMisesPlastic::correct
(
    surfaceSymmTensorField& sigma
)
{
    if (mesh().foundObject<surfaceTensorField>("grad(DD)f"))
    {
        // Lookup gradient of displacement increment
        const surfaceTensorField& gradDD =
            mesh().lookupObject<surfaceTensorField>("grad(DD)f");

        // Update the relative deformation gradient
        relFf() = I + gradDD.T();
    }
    else
    {
        // Lookup gradient of displacement
        const surfaceTensorField& gradD =
            mesh().lookupObject<surfaceTensorField>("grad(D)f");

        // Calculate the gradient of the displacement increment
        const surfaceTensorField gradDD = gradD - gradD.oldTime();

        // Update the relative deformation gradient
        relFf() = I + gradDD.T();
    }

    // Store field to calculate residual
    kappaf_.storePrevIter();

    // Take references to the internal fields for efficiency
    const tensorField& relFI = relFf().internalField();
    const symmTensorField& bEbarOldI = bEbarf_.oldTime().internalField();
    const scalarField& JOldI = Jf().oldTime().internalField();
    const scalarField& kappaOldI = kappaf_.oldTime().internalField();
    symmTensorField& bEbarI = bEbarf_.internalField();
    scalarField& JI = Jf().internalField();
    scalarField& kappaI = kappaf_.internalField();
    symmTensorField& sigmaI = sigma.internalField();

    // Calculate stress for each internal face
    forAll(sigmaI, faceI)
    {
        calculateStress
        (
            relFI[faceI],
            bEbarOldI[faceI],
            JOldI[faceI],
            kappaOldI[faceI],
            bEbarI[faceI],
            JI[faceI],
            kappaI[faceI],
            sigmaI[faceI]
        );
    }

    // Calculate stress for each boundary face
    forAll(sigma.boundaryField(), patchI)
    {
        if (sigma.boundaryField()[patchI].type() != "empty")
        {
            // Take references to the boundary patch fields for efficiency
            const tensorField& relFP = relFf().boundaryField()[patchI];
            const symmTensorField& bEbarOldP =
                bEbarf_.oldTime().boundaryField()[patchI];
            const scalarField& JOldP = Jf().oldTime().boundaryField()[patchI];
            const scalarField& kappaOldP =
                kappaf_.oldTime().boundaryField()[patchI];
            symmTensorField& bEbarP = bEbarf_.boundaryField()[patchI];
            scalarField& JP = Jf().boundaryField()[patchI];
            scalarField& kappaP = kappaf_.boundaryField()[patchI];
            symmTensorField& sigmaP = sigma.boundaryField()[patchI];

            forAll(sigmaP, faceI)
            {
                calculateStress
                (
                    relFP[faceI],
                    bEbarOldP[faceI],
                    JOldP[faceI],
                    kappaOldP[faceI],
                    bEbarP[faceI],
                    JP[faceI],
                    kappaP[faceI],
                    sigmaP[faceI]
                );
            }
        }
    }

    // Correct coupled boundaries
    bEbarf_.correctBoundaryConditions();
    Jf().correctBoundaryConditions();
    kappaf_.correctBoundaryConditions();
    sigma.correctBoundaryConditions();
}


Foam::scalar Foam::neoHookeanElasticSmoothMisesPlastic::residual()
{
    // Calculate residual based on change in plastic strain increment
    if
    (
        mesh().db().lookupObject<fvMesh>
        (
            baseMeshRegionName()
        ).foundObject<surfaceTensorField>("Ff")
    )
    {
        return
            gMax
            (
                mag
                (
                    kappaf_.internalField()
                  - kappaf_.prevIter().internalField()
                )
            )/gMax(SMALL + mag(kappaf_.prevIter().internalField()));
    }
    else
    {
        return
            gMax
            (
                mag
                (
                    kappa_.internalField()
                  - kappa_.prevIter().internalField()
                )
            )/gMax(SMALL + mag(kappa_.prevIter().internalField()));
    }
}


void Foam::neoHookeanElasticSmoothMisesPlastic::updateTotalFields()
{
    Info<< nl << "Updating total accumulated fields" << endl;

    // Reset active yield field to zero
    activeYield_ = 0.0;

    // Count cells actively yielding

    int numCellsYielding = 0;

    // Calculate increment of kappa to check which cells are yielding
    const volScalarField deltaKappa = mag(kappa_ - kappa_.oldTime());
    const volScalarField deltaKappaAvf =
        fvc::average(mag(kappaf_ - kappaf_.oldTime()));

    const scalarField& deltaKappaI = deltaKappa.internalField();
    const scalarField& deltaKappaAvfI = deltaKappaAvf.internalField();
    scalarField& activeYieldI = activeYield_.internalField();

    forAll(activeYieldI, cellI)
    {
        if (deltaKappaI[cellI] > SMALL || deltaKappaAvfI[cellI] > SMALL)
        {
            activeYieldI[cellI] = 1.0;
            numCellsYielding++;
        }
    }

    reduce(numCellsYielding, sumOp<int>());

    forAll(activeYield_.boundaryField(), patchI)
    {
        if
        (
            !activeYield_.boundaryField()[patchI].coupled()
         && activeYield_.boundaryField()[patchI].type() != "empty"
        )
        {
            const scalarField& deltaKappaP = deltaKappa.boundaryField()[patchI];
            const scalarField& deltaKappaAvfP =
                deltaKappaAvf.boundaryField()[patchI];
            scalarField& activeYieldP = activeYield_.boundaryField()[patchI];

            forAll(activeYieldP, faceI)
            {
                if (deltaKappaP[faceI] > SMALL || deltaKappaAvfP[faceI] > SMALL)
                {
                    activeYieldP[faceI] = 1.0;
                }
            }
        }
    }

    activeYield_.correctBoundaryConditions();

    Info<< "    " << numCellsYielding << " cells are actively yielding"
        << nl << endl;


    if (mesh().time().outputTime())
    {
        Info<< "Writing det(bEbar)" << endl;

        volScalarField detBEbarMinus1
        (
            "detBEbarMinus1",
            det(bEbar_) - 1.0
        );

        detBEbarMinus1.write();
    }
}


// ************************************************************************* //
