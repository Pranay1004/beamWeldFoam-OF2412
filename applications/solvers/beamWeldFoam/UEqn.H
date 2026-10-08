    MRF.correctBoundaryVelocity(U);

    const dimensionedScalar DarcyConstantlarge
    (
        "DarcyConstantlarge", dimensionSet(1,-3,-1,0,0), scalar(1.0e6)
    ); //1e10
    const dimensionedScalar DarcyConstantsmall
    (
        "DarcyConstantsmall",
        dimensionSet(0,0,0,0,0),
        scalar(1.0e-12)
    ); //1e-6

    const dimensionedScalar TSmooth("TSmooth",dimensionSet(0, 0, 0, 1, 0),50.0);

    DC =
        (
            DarcyConstantlarge*Foam::pow((1.0 - epsilon1mask), 2)
           /(Foam::pow(epsilon1mask, 3) + DarcyConstantsmall)
        );

    const volVectorField gradAlpha(fvc::grad(alpha1));

    const volVectorField nHatM(gradAlpha/(mag(gradAlpha) + deltaN));
    Marangoni = Marangoni_Constant*(gradT-nHatM*(nHatM & gradT))*mag(gradAlpha);

    if(damperSwitch)
    {
        damper = 2.0 * rho /(rho1+rho2); //To account for high density difference
    }

    // Tpos: see TEqn.H - keeps the recoil-pressure exponential defined
    // if an explicit overshoot ever drives a cell to T <= 0 (no change
    // for physical cells).
    const volScalarField TposU
    (
        max
        (
            Temperature,
            dimensionedScalar("TposU", Temperature.dimensions(), SMALL)
        )
    );

    // recoilScale (MELTING dict, default 1 = upstream behaviour): the
    // equilibrium vapour recoil 0.54*p0 overestimates digging for an arc
    // welding pool sitting in a ~1 atm plasma, so only a fraction is
    // applied as mechanical load. The evaporative heat sink Qv in TEqn.H
    // is intentionally left untouched, so the temperature field (and the
    // recoil evaluation point on the p_sat curve) does not move.
    scalar recoilScale(meltingDict.lookupOrDefault<scalar>("recoilScale", 1.0));

    pVap = min
    (
        max((Temperature - (Tvap - (TSmooth/2.0)))/TSmooth, 0.0),
        1.0
    )*0.54*p0*Foam::exp
    (
        LatentHeatVap*Mm*((TposU-Tvap)/(R*TposU*Tvap))
    )*recoilScale;

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(/*fvc::interpolate(epsilon1mask)*/rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
      + fvm::Sp(DC, U)
      - Marangoni*damper
     ==
        fvOptions(rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture.surfaceTensionForce()*fvc::interpolate(damper)
                  + fvc::interpolate(pVap)*fvc::snGrad(alpha1)
                   *fvc::interpolate(damper)
                  - ghf*fvc::snGrad(rho*rhok)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
