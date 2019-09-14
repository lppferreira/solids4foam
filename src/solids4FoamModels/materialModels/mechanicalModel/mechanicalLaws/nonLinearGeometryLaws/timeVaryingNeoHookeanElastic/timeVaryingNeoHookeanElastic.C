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

#include "timeVaryingNeoHookeanElastic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(timeVaryingNeoHookeanElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, timeVaryingNeoHookeanElastic, nonLinGeomMechLaw
    );
}


// * * * * * * * * * * *  Private Member Funtcions * * * * * * * * * * * * * //

void Foam::timeVaryingNeoHookeanElastic::makeF()
{
    if (FPtr_)
    {
        FatalErrorIn("void Foam::timeVaryingNeoHookeanElastic::makeF()")
            << "pointer already set" << abort(FatalError);
    }

    FPtr_ =
        new volTensorField
        (
            IOobject
            (
                "F",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedTensor("I", dimless, I)
        );
}


Foam::volTensorField& Foam::timeVaryingNeoHookeanElastic::F()
{
    if (!FPtr_)
    {
        makeF();
    }

    return *FPtr_;
}


void Foam::timeVaryingNeoHookeanElastic::makeFf()
{
    if (FfPtr_)
    {
        FatalErrorIn("void Foam::timeVaryingNeoHookeanElastic::makeFf()")
            << "pointer already set" << abort(FatalError);
    }

    FfPtr_ =
        new surfaceTensorField
        (
            IOobject
            (
                "Ff",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedTensor("I", dimless, I)
        );
}


Foam::surfaceTensorField& Foam::timeVaryingNeoHookeanElastic::Ff()
{
    if (!FfPtr_)
    {
        makeFf();
    }

    return *FfPtr_;
}


void Foam::timeVaryingNeoHookeanElastic::correctTimeVaryingElasticParameters()
{
    // Update elastic parameters if the time varying formulation is used
    if (timeVaryingElasticity_)
    {
        // Calculate time varying Young's modulus
        if (mesh().time().value() < transitionPeriod_ && transitionPeriod_ > SMALL)
        {
            scalar current_time = mesh().time().value();
            E_ = (desiredE_ - initialE_) * current_time / transitionPeriod_ + initialE_;
        }
        else
        {
            E_ = desiredE_;
            timeVaryingElasticity_ = false;
        }

        // Set the shear modulus
        mu_ = E_/(2.0*(1.0 + nu_));

        // Set the bulk modulus
        K_ =
        (
            planeStress()
          ? (nu_*E_/((1.0 + nu_)*(1.0 - nu_))) + (2.0/3.0)*mu_
          : (nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))) + (2.0/3.0)*mu_
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::timeVaryingNeoHookeanElastic::timeVaryingNeoHookeanElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    rho_(dict.lookup("rho")),
    timeVaryingElasticity_(false),
    transitionPeriod_(0.0),
    initialE_("initialE", dimPressure, 0.0),
    desiredE_("desiredE", dimPressure, 0.0),
    E_("E", dimPressure, 0.0),
    nu_("nu", dimless, 0.0),
    mu_("mu", dimPressure, 0.0),
    K_("K", dimPressure, 0.0),
    FPtr_(NULL),
    FfPtr_(NULL)
{
    // Read elastic parameters
    // The user can specify initial and desired E and transition period or E and nu or mu and K
    if
    (
        dict.found("initialE") && dict.found("desiredE")
     && dict.found("transitionPeriod") && dict.found("nu")
    )
    {
        // Read the Young's modulus and transition period
        timeVaryingElasticity_ = true;
        transitionPeriod_ = readScalar(dict.lookup("transitionPeriod"));
        initialE_ = dimensionedScalar(dict.lookup("initialE"));
        desiredE_ = dimensionedScalar(dict.lookup("desiredE"));

        // Calculate time varying Young's modulus
        if (mesh.time().value() < transitionPeriod_ && transitionPeriod_ > SMALL)
        {
            scalar current_time = mesh.time().value();
            E_ = (desiredE_ - initialE_) * current_time / transitionPeriod_ + initialE_;
        }
        else
        {
            E_ = desiredE_;
        }

        // Read the Poisson's ratio
        nu_ = dimensionedScalar(dict.lookup("nu"));
    }
    else if (dict.found("E") && dict.found("nu"))
    {
        // Read the Young's modulus
        E_ = dimensionedScalar(dict.lookup("E"));

        // Read the Poisson's ratio
        nu_ = dimensionedScalar(dict.lookup("nu"));
    }
    else
    {
        FatalErrorIn
        (
            "timeVaryingNeoHookeanElastic::timeVaryingNeoHookeanElastic::()"
        )   << "Either initial and desired E and transition period or E and nu "
            << "should be specified" << abort(FatalError);
    }

    // Set the shear modulus
    mu_ = E_/(2.0*(1.0 + nu_));

    // Set the bulk modulus
    K_ =
    (
        planeStress()
      ? (nu_*E_/((1.0 + nu_)*(1.0 - nu_))) + (2.0/3.0)*mu_
      : (nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))) + (2.0/3.0)*mu_
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeVaryingNeoHookeanElastic::~timeVaryingNeoHookeanElastic()
{
    deleteDemandDrivenData(FPtr_);
    deleteDemandDrivenData(FfPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::timeVaryingNeoHookeanElastic::rho() const
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


Foam::tmp<Foam::volScalarField> Foam::timeVaryingNeoHookeanElastic::impK() const
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
            (4.0/3.0)*mu_ + K_ // == 2*mu + lambda
        )
    );
}


void Foam::timeVaryingNeoHookeanElastic::correct(volSymmTensorField& sigma)
{
    // Update elastic parameters
    correctTimeVaryingElasticParameters();

    // Check if the mathematical model is in total or updated Lagrangian form
    if (nonLinGeom() == nonLinearGeometry::UPDATED_LAGRANGIAN)
    {
        if (!incremental())
        {
            FatalErrorIn(type() + "::correct(volSymmTensorField& sigma)")
                << "Not implemented for non-incremental updated Lagrangian"
                << abort(FatalError);
        }

        // Lookup gradient of displacement increment
        const volTensorField& gradDD =
            mesh().lookupObject<volTensorField>("grad(DD)");

        // Calculate the relative deformation gradient
        const volTensorField relF = I + gradDD.T();

        // Update the total deformation gradient
        F() = relF & F().oldTime();

        if (enforceLinear())
        {
            WarningIn
            (
                "void Foam::timeVaryingNeoHookeanElastic::"
                "correct(volSymmTensorField& sigma)"
            )   << "Material linearity enforced for stability!" << endl;

            // Calculate stress using Hooke's law
            sigma =
                sigma.oldTime()
              + 2.0*mu_*symm(gradDD) + (K_ - (2.0/3.0)*mu_)*tr(gradDD)*I;

            return;
        }
    }
    else if (nonLinGeom() == nonLinearGeometry::TOTAL_LAGRANGIAN)
    {
        if (incremental())
        {
            // Lookup gradient of displacement increment
            const volTensorField& gradDD =
                mesh().lookupObject<volTensorField>("grad(DD)");

            // Update the total deformation gradient
            // Note: grad is wrt reference configuration
            F() = F().oldTime() + gradDD.T();

            // Update the relative deformation gradient: not needed
            //relF() = F() & inv(F().oldTime());

            if (enforceLinear())
            {
                WarningIn
                (
                    "void Foam::timeVaryingNeoHookeanElastic::"
                    "correct(volSymmTensorField& sigma)"
                )   << "Material linearity enforced for stability!" << endl;

                // Calculate stress using Hooke's law
                sigma =
                    sigma.oldTime()
                  + 2.0*mu_*dev(symm(gradDD)) + K_*tr(gradDD)*I;

                return;
            }
        }
        else
        {
            // Lookup gradient of displacement
            const volTensorField& gradD =
                mesh().lookupObject<volTensorField>("grad(D)");

            // Update the total deformation gradient
            F() = I + gradD.T();

            // Update the relative deformation gradient: not needed
            //relF() = F() & inv(F().oldTime());

            if (enforceLinear())
            {
                WarningIn
                (
                    "void Foam::timeVaryingNeoHookeanElastic::"
                    "correct(volSymmTensorField& sigma)"
                )   << "Material linearity enforced for stability!" << endl;

                // Calculate stress using Hooke's law
                sigma = 2.0*mu_*dev(symm(gradD)) + K_*tr(gradD)*I;

                return;
            }
        }
    }
    else
    {
        FatalErrorIn
        (
            "void Foam::timeVaryingNeoHookeanElastic::"
            "correct(volSymmTensorField& sigma)"
        )   << "Unknown nonLinGeom type: " << nonLinGeom() << abort(FatalError);
    }

    // Calculate the Jacobian of the deformation gradient
    const volScalarField J = det(F());

    // Calculate the volume preserving left Cauchy Green strain
    const volSymmTensorField bEbar = pow(J, -2.0/3.0)*symm(F() & F().T());

    // Calculate the deviatoric stress
    const volSymmTensorField s = mu_*dev(bEbar);

    // Calculate the Cauchy stress
    sigma = (1.0/J)*(0.5*K_*(pow(J, 2) - 1)*I + s);
}


void Foam::timeVaryingNeoHookeanElastic::correct(surfaceSymmTensorField& sigma)
{
    // Update elastic parameters
    correctTimeVaryingElasticParameters();

    // Check if the mathematical model is in total or updated Lagrangian form
    if (nonLinGeom() == nonLinearGeometry::UPDATED_LAGRANGIAN)
    {
        if (!incremental())
        {
            FatalErrorIn(type() + "::correct(surfaceSymmTensorField& sigma)")
                << "Not implemented for non-incremental updated Lagrangian"
                << abort(FatalError);
        }

        // Lookup gradient of displacement increment
        const surfaceTensorField& gradDD =
            mesh().lookupObject<surfaceTensorField>("grad(DD)f");

        // Update the relative deformation gradient: not needed
        const surfaceTensorField relF = I + gradDD.T();

        // Update the total deformation gradient
        Ff() = relF & Ff().oldTime();

        if (enforceLinear())
        {
            WarningIn
            (
                "void Foam::timeVaryingNeoHookeanElastic::"
                "correct(surfaceSymmTensorField& sigma)"
            )   << "Material linearity enforced for stability!" << endl;

            // Calculate stress using Hooke's law
            sigma =
                sigma.oldTime() + 2.0*mu_*dev(symm(gradDD)) + K_*tr(gradDD)*I;

            return;
        }
    }
    else if (nonLinGeom() == nonLinearGeometry::TOTAL_LAGRANGIAN)
    {
        if (incremental())
        {
            // Lookup gradient of displacement increment
            const surfaceTensorField& gradDD =
                mesh().lookupObject<surfaceTensorField>("grad(DD)f");

            // Update the total deformation gradient
            // Note: grad is wrt reference configuration
            Ff() = Ff().oldTime() + gradDD.T();

            // Update the relative deformation gradient: not needed
            //relFf() = Ff() & inv(Ff().oldTime());

            if (enforceLinear())
            {
                WarningIn
                (
                    "void Foam::timeVaryingNeoHookeanElastic::"
                    "correct(surfaceSymmTensorField& sigma)"
                )   << "Material linearity enforced for stability!" << endl;

                // Calculate stress using Hooke's law
                sigma =
                    sigma.oldTime()
                  + 2.0*mu_*dev(symm(gradDD)) + K_*tr(gradDD)*I;

                return;
            }
        }
        else
        {
            // Lookup gradient of displacement
            const surfaceTensorField& gradD =
                mesh().lookupObject<surfaceTensorField>("grad(D)f");

            // Update the total deformation gradient
            Ff() = I + gradD.T();

            // Update the relative deformation gradient: not needed
            //relF() = F() & inv(F().oldTime());

            if (enforceLinear())
            {
                WarningIn
                (
                    "void Foam::timeVaryingNeoHookeanElastic::"
                    "correct(surfaceSymmTensorField& sigma)"
                )   << "Material linearity enforced for stability!" << endl;

                // Calculate stress using Hooke's law
                sigma = 2.0*mu_*dev(symm(gradD)) + K_*tr(gradD)*I;

                return;
            }
        }
    }
    else
    {
        FatalErrorIn
        (
            "void Foam::timeVaryingNeoHookeanElastic::"
            "correct(surfaceSymmTensorField& sigma)"
        )   << "Unknown nonLinGeom type: " << nonLinGeom() << abort(FatalError);
    }

    // Calculate the Jacobian of the deformation gradient
    const surfaceScalarField J = det(Ff());

    // Calculate left Cauchy Green strain tensor with volumetric term removed
    const surfaceSymmTensorField bEbar = pow(J, -2.0/3.0)*symm(Ff() & Ff().T());

    // Calculate deviatoric stress
    const surfaceSymmTensorField s = mu_*dev(bEbar);

    // Calculate the Cauchy stress
    sigma = (1.0/J)*(0.5*K_*(pow(J, 2) - 1)*I + s);
}


// ************************************************************************* //
