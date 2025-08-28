MainPulsePropagation_04_HHG_Ar1plus_808_10nm_01.m is the source code.

Library functions in the branch
analysis directory: 1 function (main_pulse_propagation)

theory directory: 12 functions
(CapillaryAttenuationCoefficient, ElectricField_d, GroupVelocityDispersion,
 HHG_DipolePhase, InverseBremsstrahlungCoefficient, IonizationTime_long,
 IonizationTime_short, RecombinationTime, StaticIonizationRate,
 ThomsonScatteringCoefficient, TunnelingIonizationRate_Linear,
 TunnelingIonizationRate_Linear2)

### Runtime-limited curve fitting

The `main_pulse_propagation` function and the source script now use
`timed_lsqcurvefit`, a wrapper around `lsqcurvefit` configured with
`optimoptions` limits. Optional parameters control the fitting process:

- `maxIterations` – maximum solver iterations
- `maxFunctionEvaluations` – maximum function evaluations
- `maxTime` – timeout in seconds (default 120)
- `downsampleFactor` – use every N-th data point

These parameters allow data to be pre-filtered and fitting to terminate
automatically when limits are exceeded.
