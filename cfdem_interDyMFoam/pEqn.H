{
    volScalarField rUA = 1.0/UEqn.A();
    surfaceScalarField rUAf("(1|A(U))", fvc::interpolate(rUA));
    volScalarField rUAvoidfraction("(voidfraction2|A(U))",rUA*voidfraction);
    surfaceScalarField rUAfvoidfraction("(voidfraction2|A(U)F)", fvc::interpolate(rUAvoidfraction));

    U = rUA*UEqn.H();

    surfaceScalarField phiU("phiU", (fvc::interpolate(U) & mesh.Sf()));

    if (pd.needReference())
    {
        adjustPhi(phi, U, pd);
    }

    #ifdef version23
    phi = phiU + 
    (
        fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
      - ghf*fvc::snGrad(rho)
    )*rUAf*mesh.magSf();
    //+ rUAfvoidfraction*fvc::ddtCorr(U, phiByVoidfraction);
    #else
    phi = phiU + 
    (
        fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
      - ghf*fvc::snGrad(rho)
    )*rUAf*mesh.magSf();
    //+ fvc::ddtPhiCorr(rUAvoidfraction, U, phiByVoidfraction);
    #endif
    surfaceScalarField phiS(fvc::interpolate(Us) & mesh.Sf());
    phi += rUAf*(fvc::interpolate(Ksl) * phiS);


    
    if (modelType=="A")
     rUAvoidfraction = volScalarField("(voidfraction2|A(U))",rUA*voidfraction*voidfraction);

     // Update the fixedFluxPressure BCs to ensure flux consistency
    #include "fixedFluxPressureHandling.H"


    //while (pimple.correctNonOrthogonal())
    #if defined(version30)
        while (piso.correctNonOrthogonal())
    #else
        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    #endif
    {
        fvScalarMatrix pdEqn
        (
            fvm::laplacian(rUAvoidfraction, pd) == fvc::div(voidfractionf*phi) + particleCloud.ddtVoidfraction()
        );

        pdEqn.setReference(pdRefCell, pdRefValue);

	
	#if defined(version30)   
            pdEqn.solve
	    (
		mesh.solver(pd.select(piso.finalInnerIter()))	
	    );
            if (piso.finalNonOrthogonalIter())
            {
                phiByVoidfraction = phi - pdEqn.flux()/voidfractionf;
            }
	#else
            if( corr == nCorr-1 && nonOrth == nNonOrthCorr )
                #if defined(versionExt32)
                    pdEqn.solve(mesh.solutionDict().solver(pd.name() + "Final")) ; //"pdFinal"));
                #else
                    pdEqn.solve(mesh.solver(pd.name() + "Final" )); //"pdFinal"));
		    
                #endif
            else
		pdEqn.solve(mesh.solutionDict().solver(pd.name()));
                //pdEqn.solve();

            if (nonOrth == nNonOrthCorr)
            {
                phiByVoidfraction = phi - pdEqn.flux()/voidfractionf;
            }
        #endif

      } // end non-orthogonal corrector loop
    

    phi = voidfractionf*phiByVoidfraction; 
    #include "continuityErrorPhiPU.H"


    if (modelType=="B" || modelType=="Bfull")
        U -= rUA *fvc::grad(pd) - Ksl*Us*rUA;   // removed rho

        //U += rUA *fvc::reconstruct((phi-phiU)/rUAf);    //Compare both results of these two
    else
        U -= voidfraction*rUA*fvc::grad(pd) - Ksl*Us*rUA;
	//U += voidfraction*rUA *fvc::reconstruct((phi-phiU)/rUAf); //Compare both results of these two

    //Info << "what's up last pressure?" << endl;
    U.correctBoundaryConditions();

    fvOptions.correct(U);
     

    //#include "continuityErrs.H"

}
