   surfaceScalarField muEff
    (
        "muEff",
        mixture.muf()
      + fvc::interpolate(voidfractionrho*turbulence->nut())
    );


    MRF.correctBoundaryVelocity(U);
    fvVectorMatrix UEqn
    (
        rho*(fvm::ddt(voidfraction, U)) - fvm::Sp(fvc::ddt(voidfractionrho),U)
      + fvm::div(rhoPhi, U) - fvm::Sp(fvc::div(rhoPhi),U)
      - fvm::laplacian(muEff, U) 
      //+ particleCloud.divVoidfractionTau(U, voidfractionrho) 
       ==
      - fvm::Sp(Ksl,U)
      + fvOptions(rho, U) 
    );	

    UEqn.relax();
    fvOptions.constrain(UEqn);

    #if defined(version30)
	if (pimple.momentumPredictor())
    #else
        if (momentumPredictor)
    #endif
    {
	if (modelType=="B" || modelType=="Bfull")
	{
	    solve
            (
	        UEqn == 
	        fvc::reconstruct
                (
                   (
			  mixture.surfaceTensionForce()
                  	- ghf*fvc::snGrad(rho)
                  	- fvc::snGrad(pd)
                    )*mesh.magSf()
                ) + Ksl*Us
	    );
            fvOptions.correct(U);
	    
	}

	else
	{
	    solve
	    (
		UEqn == 

		voidfraction*
		fvc::reconstruct
                (
                    (
			  mixture.surfaceTensionForce()
                  	- ghf*fvc::snGrad(voidfractionrho)
                  	- fvc::snGrad(pd)
                    )*mesh.magSf()

		) + Ksl*Us
	    );
   	    fvOptions.correct(U);
	}

    }
