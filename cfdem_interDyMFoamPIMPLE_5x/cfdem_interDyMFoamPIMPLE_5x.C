/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright (C) 1991-2009 OpenCFD Ltd.
                                Copyright (C) 2009-2012 JKU, Linz
                                Copyright (C) 2012-     DCS Computing GmbH,Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling.  If not, see <http://www.gnu.org/licenses/>.

Application
    cfdemSolverPiso

Description
    Transient solver for incompressible flow.
    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.
    The code is an evolution of the solver pisoFoam in OpenFOAM(R) 1.6,
    where additional functionality for CFD-DEM coupling is added.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
//#include "MULES.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "subCycle.H"
#include "CrankNicolsonDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "pimpleControl.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "interfaceProperties.H"
#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "backwardDdtScheme.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "SLTSDdtScheme.H"
#include <assert.h> 


#include "OFversion.H"
#if defined(version30)
    #include "turbulentTransportModel.H"
#else
    #include "turbulenceModel.H"
#endif
#if defined(versionv1606plus) || defined(version40)
    #include "fvOptions.H"
#else
    #include "fvIOoptionList.H"
#endif
#include "fixedFluxPressureFvPatchScalarField.H"
#include "cfdemCloud.H"

#if defined(anisotropicRotation)
    #include "cfdemCloudRotation.H"
#endif
#if defined(superquadrics_flag)
    #include "cfdemCloudRotationSuperquadric.H"
#endif
#include "implicitCouple.H"
#include "clockModel.H"
#include "smoothingModel.H"
#include "forceModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"    

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"

    #if defined(version30)
	//pimpleControl pimple(mesh);
        #include "createTimeControls.H"
    #endif
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createAlphaFluxes.H"
    #include "createFvOptions.H"
    #include "correctPhi.H"

    // create cfdemCloud
    //#include "readGravitationalAcceleration.H"
    #include "checkImCoupleM.H"
    #if defined(anisotropicRotation)
        cfdemCloudRotation particleCloud(mesh);
    #elif defined(superquadrics_flag)
        cfdemCloudRotationSuperquadric particleCloud(mesh);
    #else
        cfdemCloud particleCloud(mesh);
    #endif
    #include "checkModelType.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

  	//if (LTS)
        {
           // #include "setRDeltaT.H"
        }
        //else
        {
        #if defined(version30)
            #include "readTimeControls.H"
	    #include "alphaCourantNo.H"
 	    #include "CourantNo.H"
            #include "setDeltaT.H"
        #else
            #include "CourantNo.H"
        #endif
	}

        // do particle stuff
        particleCloud.clockM().start(1,"Global");
        particleCloud.clockM().start(2,"Coupling");
        bool hasEvolved = particleCloud.evolve(voidfraction,Us,U);

        if(hasEvolved)
        {
            particleCloud.smoothingM().smoothenAbsolutField(particleCloud.forceM(0).impParticleForces());
        }
    
        Ksl = particleCloud.momCoupleM(particleCloud.registryM().getProperty("implicitCouple_index")).impMomSource();
        Ksl.correctBoundaryConditions();

        surfaceScalarField voidfractionf = fvc::interpolate(voidfraction);
        phi = voidfractionf*phiByVoidfraction;

        //Force Checks
        #include "forceCheckIm.H"

        #include "solverDebugInfo.H"
        particleCloud.clockM().stop("Coupling");

        particleCloud.clockM().start(26,"Flow");

        int nCorr=(readInt(pimple.dict().lookup("nCorrectors")));
    
        int nOuterCorr =
        readInt(pimple.dict().lookup("nOuterCorrectors"));
        int nNonOrthCorr =
        readInt(pimple.dict().lookup("nNonOrthogonalCorrectors"));

        if(particleCloud.solveFlow())
        {
            // Pressure-velocity PISO corrector
            //{
		//fvc::makeRelative(phi, U);
                // Momentum predictor
 	//for (int oCorr = 0; oCorr < nOuterCorr; oCorr++)
              // --- Pressure-velocity PIMPLE corrector loop
         while (pimple.loop())
         {
       	    #include "alphaControls.H"
	    # include "alphaEqnSubCycle.H"
	    
	    mixture.correct();

            # include "UEqn.H"

            // --- PISO loop
            //#if defined(version30)
            while (pimple.correct())
	    {
	        #include "pEqn.H"
	    }
   	 }

	    mixture.correct();
            turbulence->correct();
	}
	
        else
        {
            Info << "skipping flow solution." << endl;
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        scalar maxSimTime     ( readScalar(runTime.controlDict().lookup("MaxExecution_Time" )) );
	scalar maxAlpha1      ( readScalar(runTime.controlDict().lookup("MAX_Alpha"         )) );
	scalar minAlpha1      ( readScalar(runTime.controlDict().lookup("MIN_Alpha"         )) );

 	scalar alphamax = max(alpha1).value();
	scalar alphaAvg = average(alpha1).value();

	if (runTime.elapsedClockTime() > maxSimTime || alphamax >= maxAlpha1 || alphaAvg <=minAlpha1  )
	{
	    Info    << " _____________________________________________________________________________________\n"
		    << "  ExecutionTime = " << runTime.elapsedCpuTime() << " s, "
                    << "  S1 = average(alpha) = " << alphaAvg << "  \n"
		    << "  MaxExecutionTime ("<<maxSimTime<<"),   "
                    << "  maxS1 ("<<maxAlpha1<<") or "
                    << "  minS1 ("<<minAlpha1 <<") reached, \n"
		    << "  Ending simulation. \n"
                    << "\n  Adjust boundary conditions or Re-set above keywords in file controlDict to continue. \n"
		    << " _____________________________________________________________________________________ \n"
                    << endl;
	    runTime.writeAndEnd();
	    //runTime.end();
            
	}

        particleCloud.clockM().stop("Flow");
        particleCloud.clockM().stop("Global");
    //}
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
