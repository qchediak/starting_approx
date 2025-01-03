# Code for PhD thesis.

These are the m-files associated with my thesis, tentatively title _Starting Approximations for Runge-Kutta Methods Applied to Ordinary Differential Equations and Differential-Algebraic Equations_.  The folders are organized as follows:

`pendulum_2stage_lobatto`: Solve the plane pendulum as a DAE using 2-stage Lobatto IIIA-IIIB.  Uses trivial starting approximations.

`pendulum_2stage_lobatto_SA`: Uses good starting approximations to solve the plane pendulum as a DAE over two steps using 2-stage Lobatto IIIA-IIIB, and plots the errors of the starting approximations vs. stepsize.

`pendulum_3stage_lobatto`: Solve the plane pendulum as a DAE using 3-stage Lobatto IIIA-IIIB.  Uses trivial starting approximations.

`pendulum_3stage_lobatto_SA`: Uses good starting approximations to solve the plane pendulum as a DAE over two steps using 3-stage Lobatto IIIA-IIIB, and plots the errors of the starting approximations vs. stepsize.

`pendulum_ODE`: Solves the pendulum problem as an ODE.  Uses trivial starting approximations.

`pendulum_ODE_SA`: Uses good starting approximations to solve the plane pendulum as an ODE over two steps and plots the errors of the starting approximations vs. stepsize.

