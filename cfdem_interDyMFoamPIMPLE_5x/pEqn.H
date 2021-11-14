{
    volScalarField rUA = 1.0/UEqn.A();
    surfaceScalarField rUAf("(1|A(U))", fvc::interpolate(rUA));
    volScalarField rUAvoidfraction("(voidfraction2|A(U))",rUA*voidfraction);
    surfaceScalarField rUAfvoidfraction("(voidfraction2|A(U)F)", fvc::interpolate(rUAvoidfraction));


    //U = rUA*UEqn.H();


    volVectorField HbyA(constrainHbyA(rUA*UEqn.H(), U, pd));

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + fvc::interpolate(rho*rUA)*fvc::ddtCorr(U, phi)
    );
    //MRF.makeRelative(phiHbyA);
    adjustPhi(phiHbyA, U, pd);

    surfaceScalarField phig
    (
        (
            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*rUAf*mesh.magSf()
    );

    phi = phiHbyA + phig;


    surfaceScalarField phiU("phiU", (fvc::interpolate(U) & mesh.Sf()));

    surfaceScalarField phiS(fvc::interpolate(Us) & mesh.Sf());
    phi += rUAf*(fvc::interpolate(Ksl) * phiS);

    
    if (modelType=="A")
     rUAvoidfraction = volScalarField("(voidfraction2|A(U))",rUA*voidfraction*voidfraction);

     // Update the fixedFluxPressure BCs to ensure flux consistency
    #include "fixedFluxPressureHandling.H"


    //while (pimple.correctNonOrthogonal())
    //#if defined(version30)
        while (pimple.correctNonOrthogonal())
    //#else
    //    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    //#endif
    {
        fvScalarMatrix pdEqn
        (
            fvm::laplacian(rUAvoidfraction, pd) == fvc::div(voidfractionf*phi) + particleCloud.ddtVoidfraction()
        );

        pdEqn.setReference(pRefCell, pRefValue);

	
	
            pdEqn.solve
	    (
		mesh.solver(pd.select(pimple.finalInnerIter()))	
	    );
            if (pimple.finalNonOrthogonalIter())
            {
                phiByVoidfraction = phi - pdEqn.flux()/voidfractionf;
		pd.relax();
		U = HbyA + rUA*fvc::reconstruct((phig - pdEqn.flux())/rUAf) - Ksl*Us*rUA;
            	//U.correctBoundaryConditions();
            	//fvOptions.correct(U);
            }

      } // end non-orthogonal corrector loop
    

    phi = voidfractionf*phiByVoidfraction; 
    #include "continuityErrorPhiPU.H"

    U.correctBoundaryConditions();
    fvOptions.correct(U);
/*

    if (modelType=="B" || modelType=="Bfull")
        //U -= rUA *fvc::grad(pd) - Ksl*Us*rUA;   // removed rho

        U += rUA *fvc::reconstruct((phig-pdEqn.flux())/rUAf) - Ksl*Us*rUA;    //Compare both results of these two
    else
        //U -= voidfraction*rUA*fvc::grad(pd) - Ksl*Us*rUA;
	U += voidfraction*rUA *fvc::reconstruct((phi-phiU)/rUAf) - Ksl*Us*rUA;;

    U.correctBoundaryConditions();

    fvOptions.correct(U);
     */

    #include "continuityErrs.H"

    p == pd + rho*gh;

    if (pd.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        pd = p - rho*gh;
    }

}
