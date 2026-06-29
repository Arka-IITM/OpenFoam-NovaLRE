#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
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

#include "createPhi.H"


Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();


Info<< "Reading field K\n" << endl;
volScalarField K
(
    IOobject
    (
        "K",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading KProperties\n" << endl;
const dictionary& KProperties = mixture.subDict("KProperties");

dimensionedScalar rho1Intercept
(
    "rho1Intercept",
    dimMass/dimVolume,
    KProperties
);

dimensionedScalar rho1Slope
(
    "rho1Slope",
    dimMass/dimVolume,
    KProperties
);

// Liquid-phase density as a linear function of K
// rho1Field = rho1Intercept + rho1Slope*K
// K=0 -> 98/2 H2O2/water (rho1Intercept)
// K=1 -> 70/30 ethanol/water (rho1Intercept + rho1Slope)
volScalarField rho1Field
(
    IOobject
    (
        "rho1Field",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rho1Intercept + rho1Slope*K
);

dimensionedScalar mu1Eth
(
    "mu1Eth",
    (dimMass/dimLength/dimTime),
    KProperties
);

dimensionedScalar mu1H2O2
(
    "mu1H2O2",
    (dimMass/dimLength/dimTime),
    KProperties
);

// Unit-viscosity scalar, used only to non-dimensionalise mu1Eth/mu1H2O2
// before taking their log (Foam::log requires a dimensionless argument),
// then to restore Pa.s units after exp().
dimensionedScalar unitDynVisc("unitDynVisc", (dimMass/dimLength/dimTime), 1.0);

// log(mu1Eth) and log(mu1H2O2), each now dimensionless scalars
dimensionedScalar logMu1Eth(log(mu1Eth/unitDynVisc));
dimensionedScalar logMu1H2O2(log(mu1H2O2/unitDynVisc));

// Logarithmic (Grunberg-Nissan-type) DYNAMIC viscosity mixing rule,
// mass-fraction weighted by K (matches the K convention used for
// rho1Field and the species mass fractions throughout this case):
//   mu1Field = mu1Eth^K * mu1H2O2^(1-K)
//            = exp( K*log(mu1Eth) + (1-K)*log(mu1H2O2) )
// At K=1 -> mu1Eth (70% ethanol solution), at K=0 -> mu1H2O2 (98% H2O2 solution)
//
// NOTE: this produces DYNAMIC viscosity (Pa.s). OpenFOAM's "nu" transport
// model expects KINEMATIC viscosity (m^2/s), so mu1Field is divided by
// rho1Field below to get nu1Field. See Gmehling et al. - this class of
// simple log-mixing rule is known to be order-of-magnitude accurate at
// best; treat nu1Field as an approximation, not a validated correlation.
volScalarField logMu1Field
(
    IOobject
    (
        "logMu1Field",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    K*logMu1Eth + (1.0 - K)*logMu1H2O2
);

volScalarField mu1Field
(
    IOobject
    (
        "mu1Field",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    exp(logMu1Field)*unitDynVisc
);

// Kinematic viscosity = dynamic viscosity / density
volScalarField nu1Field
(
    IOobject
    (
        "nu1Field",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mu1Field/rho1Field
);

dimensionedScalar sigmaCubic
(
    "sigmaCubic",
    dimMass/dimTime/dimTime,
    KProperties
);

dimensionedScalar sigmaQuad
(
    "sigmaQuad",
    dimMass/dimTime/dimTime,
    KProperties
);

dimensionedScalar sigmaLin
(
    "sigmaLin",
    dimMass/dimTime/dimTime,
    KProperties
);

dimensionedScalar sigmaConst
(
    "sigmaConst",
    dimMass/dimTime/dimTime,
    KProperties
);

// Surface tension vs air, as a function of K, fitted to a Parachor-method
// (Sugden) bulk-mixture estimate over the full ethanol/H2O2/water ternary
// (no surface-enrichment correction - see conversation notes on why the
// full Sprow-Prausnitz/Suarez iterative model was judged not worth the
// added cost for this case). Cubic fit, max error ~0.05 mN/m vs the
// underlying Parachor calculation across K=0 to 1.
//   sigmaField(K) = sigmaCubic*K^3 + sigmaQuad*K^2 + sigmaLin*K + sigmaConst
// At K=0 -> ~98% H2O2/water, at K=1 -> 70% ethanol/water
volScalarField sigmaField
(
    IOobject
    (
        "sigmaField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    sigmaCubic*K*K*K + sigmaQuad*K*K + sigmaLin*K + sigmaConst
);


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1Field + alpha2*rho2
);
rho.oldTime();


// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);


typedef incompressibleInterPhaseTransportModel
    <
        immiscibleIncompressibleTwoPhaseMixture
    > transportModelType;

autoPtr<transportModelType> turbulence
(
    new transportModelType
    (
        rho, U, phi, rhoPhi, mixture
    )
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


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
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

#include "createMRF.H"
#include "createFvOptions.H"
#include "createPorosity.H"

isoAdvection advector(alpha1, phi, U);

dropletCloud cloud(mesh, g);

phaseCoupling coupling
(
    U,
    alpha1,
    cloud
);
