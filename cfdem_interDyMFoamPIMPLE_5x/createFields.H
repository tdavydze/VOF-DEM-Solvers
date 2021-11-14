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

   //#   include "createPhi.H"

   //#ifndef createPhi_H
    //#define createPhi_H
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
    //#endif


    Info<< "Reading transportProperties\n" << endl;

    immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

    volScalarField& alpha1(mixture.alpha1());
    volScalarField& alpha2(mixture.alpha2());

    const dimensionedScalar& rho1 = mixture.rho1();
    const dimensionedScalar& rho2 = mixture.rho2();

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
    ); 
/*
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
	//voidfraction*alpha1
	mesh
    );
*/
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

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, pd, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);


    if (pd.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        ),
        pd = p - rho*gh;
    }

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
    );

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


    mesh.setFluxRequired(pd.name());
    mesh.setFluxRequired(alpha1.name());


    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, mixture)
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
        rho1*phi   // Original
	//fvc::interpolate(rho)*phi //phi
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
	voidfraction*rho
    );

    #include "createMRF.H"

    interfaceProperties interface(alpha1, U, mixture);
