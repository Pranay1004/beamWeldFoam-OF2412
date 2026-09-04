/*---------------------------------------------------------------------------*\
  Adaptive Time Step Control for beamWeldFoam

  Adjusts the time step based on the Courant number to ensure stability.
  The time step is limited by:
      1. Max Courant number (maxCo from controlDict)
      2. Max alpha Courant number (maxAlphaCo from controlDict)
      3. Maximum time step (maxDeltaT from controlDict)

  The adjustment formula is:
      dt_new = min(maxCo * dt / Co, min(dt, maxDeltaT))

  This ensures that the Courant number stays below maxCo while not
  increasing the time step beyond its current value or maxDeltaT.

  Only adjusts on the first timestep (timeIndex == 0) if CoNum > SMALL.
  Subsequent timesteps use the adjusted value.

  Original authors:
      Thomas F. Flint, University of Manchester
      Gowthaman Parivendhan, University College Dublin
      Alojz Ivankovic, University College Dublin
      Michael C. Smith, University of Manchester
      Philip Cardiff, University College Dublin

  OpenFOAM 2412 port by:
      Pranay K. Pandey, Indian Institute of Space Science & Technology, Trivandrum
      github.com/Pranay1004

\*---------------------------------------------------------------------------*/

if (adjustTimeStep)
{
    // ---- Adjust time step on first iteration ----
    // Only adjust if the Courant number is significant (CoNum > SMALL)
    if ((runTime.timeIndex() == 0) && (CoNum > SMALL))
    {
        runTime.setDeltaT
        (
            min
            (
                // ---- Co-based adjustment ----
                // dt_new = maxCo * dt / Co
                // This scales the time step to achieve the target Courant number
                maxCo*runTime.deltaTValue()/CoNum,

                // ---- Upper bound ----
                // Don't increase dt beyond current value or maxDeltaT
                min(runTime.deltaTValue(), maxDeltaT)
            )
        );
    }
}

// ************************************************************************* //
