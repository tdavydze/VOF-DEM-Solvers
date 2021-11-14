CorrectPhi
(
    U,
    phi,
    pd,
    dimensionedScalar("rAUf", dimTime/rho.dimensions(), 1),
    geometricZeroField(),
    pimple
);

#include "continuityErrs.H"
