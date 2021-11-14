   surfaceScalarField muEff
    (
        "muEff",
        mixture.muf()
      + fvc::interpolate(rho*turbulence->nut())
    );


    MRF.correctBoundaryVelocity(U);
    fvVectorMatrix UEqn
    (
        rho*(fvm::ddt(voidfraction, U)) - fvm::Sp(fvc::ddt(voidfractionrho),U)
      + fvm::div(rhoPhi, U) - fvm::Sp(fvc::div(rhoPhi),U)  
      - fvm::laplacian(muEff, U) 
      //+ particleCloud.divVoidfractionTau(U, voidfraction)  // Still not resolved
       ==
      - fvm::Sp(Ksl,U)
	    //+ fvOptions(voidfractionrho, U)
    );	

    UEqn.relax();
    fvOptions.constrain(UEqn);

/* These guys should be corret  */
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
		//voidfraction*
	        fvc::reconstruct
                (
                   (
			 //fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1) 
			  mixture.surfaceTensionForce()
                  	- ghf*fvc::snGrad(rho)
                  	- fvc::snGrad(pd)
                    )*mesh.magSf()
                ) + Ksl*Us
	    );
            //fvOptions.correct(U);
	    
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
                   	 //fvc::interpolate(interface.sigmaK())*fvc::snGrad(volfractionAlpha)
			  mixture.surfaceTensionForce()
                  	- ghf*fvc::snGrad(voidfractionrho)
                  	- fvc::snGrad(pd)
                    )*mesh.magSf()

		) + Ksl*Us
	    );
   	    //fvOptions.correct(U);
	}

    }
