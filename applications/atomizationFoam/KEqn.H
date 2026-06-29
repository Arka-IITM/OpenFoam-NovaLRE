{
    // divPhi = discrete phase-continuity residual for the liquid phase.
    // Two mechanisms destroy the K matrix diagonal and cause SIGFPE:
    //
    //  1. fvm::ddt(alpha1, K) contributes alpha1_new/dt to the diagonal.
    //     When isoAdvection leaves alpha1 < 0 (undershoot) this is negative,
    //     and when alpha1 == 0 (pure-gas cell) it is exactly zero.
    //
    //  2. -fvm::Sp(divPhi, K): when divPhi > 0 (advancing liquid front) this
    //     subtracts from the diagonal.
    //
    // Fix strategy:
    //  a) Replace alpha1 in the ddt diagonal by max(alpha1, minAlpha) via an
    //     extra Sp term.  This guarantees the ddt diagonal >= minAlpha/dt > 0
    //     at all times.  In gas cells (alpha1 ~ 0) it acts as a K-sink driving
    //     K to zero, which is physically correct -- K has no meaning without liquid.
    //     Using minAlpha = 1e-6 gives a floor of ~4000 s^-1 at dt = 2.5e-10 s,
    //     orders of magnitude above SMALL and robust against exact cancellation
    //     between the alpha1/dt and -alpha1Neg/dt terms.
    //  b) Split divPhi: keep only its non-positive part implicit (safe) and move
    //     the positive part to the explicit RHS (bounded because K in [0,1]).

    volScalarField divPhi
    (
        "divPhi",
        fvc::ddt(alpha1) + fvc::div(phi)
    );

    const dimensionedScalar zeroPhi("zero", divPhi.dimensions(), 0.0);

    const dimensionedScalar minAlpha("minAlpha", dimless, 1e-6);
    const volScalarField alpha1Eff(max(alpha1, minAlpha));

    fvScalarMatrix KEqn
    (
        fvm::ddt(alpha1, K)
      + fvm::div(phi, K)
        // Boost ddt diagonal from alpha1/dt to max(alpha1,minAlpha)/dt.
        // (alpha1 - alpha1Eff) <= 0 always, so -Sp adds a non-negative entry.
      - fvm::Sp((alpha1 - alpha1Eff)/runTime.deltaT(), K)
        // Implicit: only the non-positive (stabilising) part of divPhi.
      - fvm::Sp(min(divPhi, zeroPhi), K)
     ==
        // Explicit: positive (destabilising) part of divPhi.
      - max(divPhi, zeroPhi) * K
      + fvOptions(K)
    );

    KEqn.relax();
    KEqn.solve();

    K.max(0.0);
    K.min(1.0);
}
