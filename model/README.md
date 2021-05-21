Directory `model/` contains the core supraglacial hydrology model code.
The files in this directory are:

# Directory contents

## ODE integrators:
 * `odeFE.m`: Function that provides simple first-order Forward Euler time stepping. The usage is identical to built-in MATLAB ode integrators (e.g. ODE45) except this function also takes a timestep argument (`dt`).
 * `odeRK.m`: Function that provides four steo Runge-Kutta integration (RK4 method). The usage is identical to `odeFE.m`.

## Helper functions
 * `get_default_params.m` returns a parameter structure containing all the model default parameters. This is also the central place that stores **all** default parameters. Defaults should be stored here, and only here.
 * `get_derivatives.m` runs the model (e.g. `rhs_sheet_and_channel.m`) but returns the time derivatives and secondary variables (e.g. flux and potentials) at a single time instead of integrating over a range of times.
 * `validate_params.m`.

## Model code
 * `rhs_sheet_and_channel.m` computes the time derivative (e.g. Right-hand side) of the model equatiions. This is the core of the numerical model. You won't usually interact directly with this function, instead using `run_model.m` and `get_derivatives.m`.
 * `unpack_state_vector.m` unpacks the state vector into individual dynamic variables. This is a low-level function you probably won't interact with.

## Wrappers
 * `run_model.m` is the function you use to run the core model. This function chooses the time integration method (based on the provided parameters) and saves the model outputs.
 * `unpack_outputs.m` and `save_model_outputs.m` provide low-level functionality for `run_model.m`.
