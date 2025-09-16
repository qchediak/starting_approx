# Starting Approximations for Runge-Kutta Methods Applied to Ordinary Differential Equations and Differential-Algebraic Equations

These are the Matlab files associated with my thesis, tentatively titled _Starting Approximations for Runge-Kutta Methods Applied to Ordinary Differential Equations and Differential-Algebraic Equations_.  These are meant to illustrate and verify the analytical results of the thesis.  This code could be modified for other problems; however, these scripts are not written as general solvers and are not ready for operational use. 

Each folder is self-contained (i.e., code in one directory does not reference code in another directory).  To run any of the experiments, navigate within a folder and run `main` in MATLAB.  To change parameters, look in the `init` or `main` files.  The code in this repository was run in MATLAB R2024b; earlier versions may have formatting issues in the plots.

The folders are organized as follows:

- `pendulum_2stage_lobatto`: Solve the plane pendulum as a DAE using 2-stage Lobatto IIIA-IIIB.  Uses trivial starting approximations.

- `pendulum_2stage_lobatto_SA`: Uses good starting approximations to solve the plane pendulum as a DAE over two steps using 2-stage Lobatto IIIA-IIIB, and plots the errors of the starting approximations vs. stepsize.

- `pendulum_3stage_lobatto`: Solve the plane pendulum as a DAE using 3-stage Lobatto IIIA-IIIB.  Uses trivial starting approximations.

- `pendulum_3stage_lobatto_SA`: Uses good starting approximations to solve the plane pendulum as a DAE over two steps using 3-stage Lobatto IIIA-IIIB, and plots the errors of the starting approximations vs. stepsize.

- `pendulum_ODE`: Solves the pendulum problem as an ODE.  Uses trivial starting approximations.

- `pendulum_ODE_SA`: Uses good starting approximations to solve the plane pendulum as an ODE over two steps and plots the errors of the starting approximations vs. stepsize.

- `two_body_ode`: Solves the two-body problem using trivial starting approximations.

- `two_body_ode_SA`: Uses good starting approximations to solve the two-body problem over two steps and plots the errors of the starting approximations vs. stepsize.

- `pendulum_PRK`: Solves the plane pendulum as a DAE using PRK methods implemented in the standard way.  Uses trivial starting approximations.

- `pendulum_PRK_SA`: Uses good starting approximations to solve the plane pendulum as a DAE over two steps using PRK methods implemented in the standard way.  Also plots the errors of the starting approximations vs. stepsize.

