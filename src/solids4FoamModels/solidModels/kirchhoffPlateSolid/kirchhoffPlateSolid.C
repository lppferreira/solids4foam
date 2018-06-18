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

#include "kirchhoffPlateSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "faCFD.H"
#include "linearElastic.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kirchhoffPlateSolid, 0);
addToRunTimeSelectionTable(physicsModel, kirchhoffPlateSolid, solid);
addToRunTimeSelectionTable(solidModel, kirchhoffPlateSolid, dictionary);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

bool kirchhoffPlateSolid::converged
(
    const int iCorr,
    const lduSolverPerformance& solverPerfM,
    const lduSolverPerformance& solverPerfw,
    const areaScalarField& M,
    const areaScalarField& w
)
{
    // We will check a number of different residuals for convergence
    bool converged = false;

    // Calculate relative residuals
    const scalar residualM =
        gMax
        (
            mag(M.internalField() - M.prevIter().internalField())
           /max
            (
                gMax(mag(M.internalField() - M.oldTime().internalField())),
                SMALL
            )
        );

    const scalar residualw =
        gMax
        (
            mag(w.internalField() - w.prevIter().internalField())
           /max
            (
                gMax(mag(w.internalField() - w.oldTime().internalField())),
                SMALL
            )
        );

    // Calculate material residual
    const scalar materialResidual = mechanical().residual();

    // If one of the residuals has converged to an order of magnitude
    // less than the tolerance then consider the solution converged
    // force at leaast 1 outer iteration and the material law must be converged
    if (iCorr > 1 && materialResidual < materialTol())
    {
        bool convergedM = false;
        bool convergedw = false;

        if
        (
            (
                solverPerfM.initialResidual() < solutionTol()
             && residualM < solutionTol()
            )
         || solverPerfM.initialResidual() < alternativeTol()
         || residualM < alternativeTol()
        )
        {
            convergedM = true;
        }

        if
        (
            (
                solverPerfw.initialResidual() < solutionTol()
             && residualw < solutionTol()
            )
         || solverPerfw.initialResidual() < alternativeTol()
         || residualw < alternativeTol()
        )
        {
            convergedw = true;
        }


        if (convergedM && convergedw)
        {
            Info<< "    The residuals have converged" << endl;
            converged = true;
        }
    }

    // Print residual information
    if (iCorr == 0)
    {
        Info<< "    Corr, res (M & w), relRes (M & w), matRes, iters (M & w)"
            << endl;
    }
    else if (iCorr % infoFrequency() == 0 || converged)
    {
        Info<< "    " << iCorr
            << ", " << solverPerfM.initialResidual()
            << ", " << solverPerfw.initialResidual()
            << ", " << residualM
            << ", " << residualw
            << ", " << materialResidual
            << ", " << solverPerfM.nIterations()
            << ", " << solverPerfw.nIterations() << endl;

        if (converged)
        {
            Info<< endl;
        }
    }
    else if (iCorr == nCorr() - 1)
    {
        maxIterReached()++;
        Warning
            << "Max iterations reached within the M-w loop" << endl;
    }

    return converged;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kirchhoffPlateSolid::kirchhoffPlateSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    aMesh_(mesh()),
    w_
    (
        IOobject
        (
            "w",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh_
    ),
    wVf_
    (
        IOobject
        (
            "wVf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimLength, 0.0)
    ),
    M_
    (
        IOobject
        (
            "M",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh_
    ),
    MVf_
    (
        IOobject
        (
            "MVf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimPressure/dimArea, 0.0)
    ),
    p_
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh_
    ),
    theta_
    (
        IOobject
        (
            "theta",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh_,
        dimensionedVector("zero", dimless, vector::zero)
    ),
    thetaVf_
    (
        IOobject
        (
            "thetaVf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", dimless, vector::zero)
    ),
    gradTheta_(fac::grad(theta_)),
    rho_("zero", dimDensity, 0.0),
    E_("zero", dimPressure, 0.0),
    nu_("zero", dimless, 0.0),
    h_(solidProperties().lookup("plateThickness")),
    bendingStiffness_("zero", dimPressure*dimVolume, 0.0),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_)
{
    const PtrList<mechanicalLaw>& mechLaws = mechanical();

    // Only the linearElastic mechanicalLaw is allow and one material
    if (mechLaws.size() != 1)
    {
        FatalErrorIn(type() + "::" + type())
            << " can currently only be used with a single material"
            << abort(FatalError);
    }
    else if (!isA<linearElastic>(mechLaws[0]))
    {
        FatalErrorIn(type() + "::" + type())
            << " can only be used with the linearElastic "
            << "mechanicalLaw" << nl
            << abort(FatalError);
    }

    // Cast the mechanical law to a linearElastic mechanicalLaw
    const linearElastic& mech = refCast<const linearElastic>(mechLaws[0]);

    // Set plate properties
    rho_ = mech.rhoScalar();
    E_ = mech.E();
    nu_ = mech.nu();
    bendingStiffness_ = E_*pow(h_, 3)/(12*(1 - pow(nu_, 2)));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool kirchhoffPlateSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    // Create volume-to surface mapping object
    volSurfaceMapping vsm(aMesh_);

    // Mesh update loop
    do
    {
        int iCorr = 0;
        lduSolverPerformance solverPerfM;
        lduSolverPerformance solverPerfw;
        blockLduMatrix::debug = 0;

        Info<< "Solving the Kirchhoff plate equation for w and M" << endl;

        // w and M equation loop
        do
        {
            // Algorithm
            // Solve M equation
            // Solve w equation
            // where
            // M is the moment sum
            // w is the transvere (out of plane) displacement
            // The M equation is:
            //     rho*h*fac::d2dt2(w) = fam::laplacian(M) + p
            // and the w equation is:
            //     fam::laplacian(D, w) + M
            // where
            // rho is the density
            // h is the plate thickness
            // p is the net transverse pressure
            // D is the bending stiffness
            // THESE ARE BOTH SCALARS: IDEA: USE BLOCK COUPLED!

            // Store fields for under-relaxation and residual calculation
            M_.storePrevIter();

            // Solve M equation
            // d2dt2 not implemented so calculate it manually using ddt
            // Also, "==" complains so we will move all terms to left
            faScalarMatrix MEqn
            (
                rho_*h_
               *(
                   fac::ddt(w_) - fac::ddt(w_.oldTime())
                )/runTime().deltaT()
              - fam::laplacian(M_) - p_
            );

            // Relax the linear system
            MEqn.relax();

            // Solve the linear system
            solverPerfM = MEqn.solve();

            // Relax the field
            M_.relax();

            // Store fields for under-relaxation and residual calculation
            w_.storePrevIter();

            // Solve w equation
            faScalarMatrix wEqn
            (
                fam::laplacian(bendingStiffness_, w_) + M_
            );

            // Relax the linear system
            wEqn.relax();

            // Solve the linear system
            solverPerfw = wEqn.solve();

            // Relax the field
            w_.relax();

            // Update the angle of rotation
            theta_ = -fac::grad(w_);

            // Update the gradient of rotation field, used for non-orthogonal
            // correction in clamped boundary conditions
            gradTheta_ = fac::grad(theta_);
        }
        while
        (
            !converged(iCorr, solverPerfM, solverPerfw, M_, w_)
         && ++iCorr < nCorr()
        );

        // Map area fields to vol fields
        vsm.mapToVolume(M_, MVf_.boundaryField());
        vsm.mapToVolume(w_, wVf_.boundaryField());
        vsm.mapToVolume(theta_, thetaVf_.boundaryField());

        // Set the z component of the displacement field (we are assuming the
        // plate normal is in the z direction)
        D() = wVf_*vector(0, 0, 1);

        // Interpolate cell displacements to vertices
        mechanical().interpolate(D(), pointD());

        // Increment of displacement
        DD() = D() - D().oldTime();

        // Increment of point displacement
        pointDD() = pointD() - pointD().oldTime();

        // Velocity
        U() = fvc::ddt(D());
    }
    while (mesh().update());

    return true;
}


tmp<vectorField> kirchhoffPlateSolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch mechanical property
    const scalarField& impK = impK_.boundaryField()[patchID];

    // Patch reciprocal implicit stiffness field
    const scalarField& rImpK = rImpK_.boundaryField()[patchID];

    // Patch gradient
    const tensorField& pGradD = gradD().boundaryField()[patchID];

    // Patch stress
    const symmTensorField& pSigma = sigma().boundaryField()[patchID];

    // Patch unit normals
    const vectorField n = patch.nf();

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - n*pressure)
              - (n & (pSigma - impK*pGradD))
            )*rImpK
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
