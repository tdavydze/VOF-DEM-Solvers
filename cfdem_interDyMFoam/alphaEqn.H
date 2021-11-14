{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phic = mag(phi/mesh.magSf());
    phic = min(interface.cAlpha()*phic, max(phic));
    surfaceScalarField phir = phic*interface.nHatf();

    forAll (mesh.cells(), cellI)
    {
        meshVol = mesh.V()[cellI];
	invMeshVol = scalar(1.0)/meshVol;  
        phase1Mass = rho1alpha1*meshVol;
    }
    Info << " Running void fraction Alpha  "   << endl;
    volfractionAlpha = (phase1Mass/rho1)*invMeshVol;


    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        surfaceScalarField phiAlpha =
            fvc::flux
            (
                phi,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, scalar(1) - alpha1, alpharScheme),
                alpha1,
                alpharScheme
            );

        MULES::explicitSolve(alpha1, phi, phiAlpha, 1, 0);
	/*surfaceScalarField phiAlpha =
            fvc::flux
            (
                phi,
                volfractionAlpha, //alpha1,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, scalar(1) - volfractionAlpha, alpharScheme),
                volfractionAlpha,  //voidfraction*alpha1,
                alpharScheme
            );

        MULES::explicitSolve(volfractionAlpha, phi, phiAlpha, 1, 0);

         
        rhoPhi = phiAlpha*(rho1 - rho2) + phi*rho2;  //mass flux 
       	//rhophi2 = phiAlpha*(rho2 - rho1) + phi*rho1;     for secondary fluid*/
    }

    Info<< "Liquid phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Max(rhoPhi) = " << max(rhoPhi).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value()
        << endl;
}
