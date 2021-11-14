    Info<< "Reading field p\n" << endl;
    
    volScalarField pd
    (
        IOobject
        (
            "pd",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    Info<< "Reading physical velocity field U" << endl;
    Info<< "Note: only if voidfraction at boundary is 1, U is superficial velocity!!!\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//================================
// Two phase flow properties
//================================
    volScalarField alpha1
    (
	IOobject
	(
	    "alpha1",
	    runTime.timeName(),
	    mesh,
	    IOobject::MUST_READ,
	    IOobject::AUTO_WRITE
	),
	mesh
    );
    
    //alpha1.boundaryField().evaluate();

   Info<< "\nReading voidfraction field voidfraction = (Vgas/Vparticle)\n" << endl;
    volScalarField voidfraction
    (
        IOobject
        (
            "voidfraction",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
     
    /*tmp<volScalarField> meshVol
    {
	return tmp<volScalarField>
	(
	    new volScalarField("MeshVolume", mesh.V())
	);
    }*/


    volScalarField meshVol
    (
	IOobject
	(
	    "meshVol",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
        ),
        mesh,

	//dimensionedScalar("small", dimLength*dimLength*dimLength, SMALL)
    	//dimensionedScalar("zero",dimVolume,0.0)
	dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), 0.0)
    );
    //meshVol.internalField() = mesh.V();

    volScalarField invMeshVol
    (
	IOobject
	(
	    "invMeshVol",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
        ),
        mesh,
	//dimensionedScalar("zero",dimVolume,0.0)
        dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), 0.0)	//1.0)
    );

    
/*
    volScalarField voidfractionAlpha
    (
	IOobject
	(
	    "voidfractionAlpha",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::AUTO_WRITE
        ),
        mesh
 	//alpha1/mesh.V()
    );
*/
   #ifndef createPhi_H
    #define createPhi_H
    Info<< "Reading/calculating face flux field phi\n" << endl;
    surfaceScalarField phi
    (
	IOobject
	(
	    "phi",
	    runTime.timeName(),
	    mesh,
	    IOobject::READ_IF_PRESENT,
	    IOobject::AUTO_WRITE
	),
	linearInterpolate(U*voidfraction) & mesh.Sf()
    );
    #endif

    Info << "Reading transportProperties\n" << endl;
    immiscibleIncompressibleTwoPhaseMixture twoPhaseProperties(U, phi);

   

    dimensionedScalar rho1(twoPhaseProperties.rho1());
    dimensionedScalar rho2(twoPhaseProperties.rho2());

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        //mesh//,
        alpha1*(rho1 + (scalar(1) - alpha1) * rho2)
    );
    rho.oldTime();

    volScalarField rho1alpha1
    (
	IOobject
	(
	    "rho1alpha1",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	rho1*alpha1
	//alpha1.internalField().types()
	//alpha1.boundaryField().types()
    ); 

    volScalarField phase1Mass
    (
	IOobject
	(
	    "phaseMass",
	    runTime.timeName(),
	    mesh,
	    IOobject::MUST_READ,
	    IOobject::NO_WRITE
	),
	mesh
	//alpha1.internalField().types()
	//alpha1.boundaryField().types()
    ); 

    volScalarField volfractionAlpha
    (
	IOobject
	(
	    "volfractionAlpha",
	    runTime.timeName(),
	    mesh,
	    IOobject::MUST_READ,
	    IOobject::NO_WRITE
	),
	mesh
    );

    uniformDimensionedVectorField g
    (
        IOobject
        (
            "g",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
 
    volScalarField gh
    (
   	"gh",
	g & mesh.C()
    );

    surfaceScalarField ghf
    (
	"ghf", 
	g & mesh.Cf()
    );

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pd+rho*gh
    );

    volTensorField gradU(fvc::grad(U));

    volVectorField gradp(fvc::grad(p));
//===============================
// particle interaction modelling
//===============================

    Info<< "\nReading momentum exchange field Ksl\n" << endl;
    volScalarField Ksl
    (
        IOobject
        (
            "Ksl",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
        //dimensionedScalar("0", dimensionSet(1, -3, -1, 0, 0), 1.0)
    );


    //Info<< "\nCreating density field rho\n" << endl;
    

    Info<< "Reading particle velocity field Us\n" << endl;
    volVectorField Us
    (
        IOobject
        (
            "Us",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


//===============================

//#   include "createPhi.H"
    

Info<< "Generating interstitial face flux field phiByVoidfraction\n" << endl;
surfaceScalarField phiByVoidfraction
(
 IOobject
 (
    "phiByVoidfraction",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
 ),
 linearInterpolate(U) & mesh.Sf()
);

    label pdRefCell = 0;
    scalar pdRefValue = 0.0;
    setRefCell(pd, mesh.solutionDict().subDict("PISO"), pdRefCell, pdRefValue);


    if (pd.needReference())
    {
        pd += dimensionedScalar
        (
            "pd",
            pd.dimensions(),
            pdRefValue - getRefCellValue(pd, pdRefCell)
        ),
        pd = p - rho*gh;
    }

    mesh.setFluxRequired(pd.name());
    mesh.setFluxRequired(alpha1.name());


    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, twoPhaseProperties)
    );

    // Mass flux
    surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rho*phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        //rho1*phi   // Original
	fvc::interpolate(rho)*phi
    );

    volScalarField voidfractionrho
    (
        IOobject
        (
            "voidfraction*rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	rho*voidfraction
    );

    interfaceProperties interface(alpha1, U, twoPhaseProperties);
