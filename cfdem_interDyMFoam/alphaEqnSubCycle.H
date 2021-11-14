label nAlphaCorr
(
    readLabel(piso.dict().lookup("nAlphaCorr"))    // was pimple
);

label nAlphaSubCycles
(
    readLabel(piso.dict().lookup("nAlphaSubCycles"))
);

if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();
    surfaceScalarField rhoPhiSum = 0.0*rhoPhi;

    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
#       include "alphaEqn.H"
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
    }
    Info << "what's up phi1 "  << endl;
    rhoPhi = rhoPhiSum;
    Info << "phi2? "  << endl;
}
else
{
#       include "alphaEqn.H"
}

interface.correct();
//rho == volfractionAlpha*rho1 + (scalar(1) - volfractionAlpha)*rho2;
rho == alpha1*rho1 + (scalar(1) - alpha1)*rho2;
