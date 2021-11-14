   surfaceScalarField muEff
    (
        "muEff",
        twoPhaseProperties.muf()
      + fvc::interpolate(rho*turbulence->nut())
    );

Info << "phi UQN? "  << endl;
    fvVectorMatrix UEqn
    (
        fvm::ddt(voidfractionrho, U) - fvm::Sp(fvc::ddt(voidfractionrho),U)
      + fvm::div(rhoPhi, U) - fvm::Sp(fvc::div(rhoPhi),U)   //voidfraction is in phi declatation.
      + particleCloud.divVoidfractionTau(U, voidfractionrho)  //should we include rho? 
       ==
      - fvm::Sp(Ksl,U)
            //+ fvOptions(U*rho)
    );	

/*
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      - fvm::laplacian(muEff, U)
      - (fvc::grad(U) & fvc::grad(muEff))
    //- fvc::div(muEff*(fvc::interpolate(dev(fvc::grad(U))) & mesh.Sf()))
    );
*/
    UEqn.relax();
    fvOptions.constrain(UEqn);

       //Info << "what's up U?" << endl;
    #if defined(version30)
	if (piso.momentumPredictor())
    #else
        if (momentumPredictor)
    #endif
    {
	if (modelType=="B" || modelType=="Bfull")
	    solve
            (
	        UEqn == 

	        fvc::reconstruct
                (
                    (
                   	 fvc::interpolate(interface.sigmaK())*fvc::snGrad(volfractionAlpha)
                  	- ghf*fvc::snGrad(rho)
                  	- fvc::snGrad(pd)
                    )*mesh.magSf()
                ) + Ksl*Us
	    );
         //Info << "what's up Uagain?" << endl;
	else
	    solve
	    (
		UEqn == 

		voidfraction*
		fvc::reconstruct
                (
                    (
                   	 fvc::interpolate(interface.sigmaK())*fvc::snGrad(volfractionAlpha)
                  	- ghf*fvc::snGrad(rho)
                  	- fvc::snGrad(pd)
                    )*mesh.magSf()

		) + Ksl*Us
	    );
	Info << "what's up Uagain?" << endl;
	fvOptions.correct(U);

    }
