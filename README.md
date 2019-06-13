# Photosynthesis Solver

Small code meant to explore the relative merits of solving the nonlinear equations of state describing photosynthesis in different ways. Need to have the `TDYCORE_DIR` environmental variable set and then run with:

```
./Photosynthesis -snes_fd -snes_monitor
```

It currently only sticks in random numbers on the unit interval for the parameters. Mostly this will solve in 3-5 interations, but occasionally it will stagnate and not solve. This is because for some parameter combinations the residual contains several poles in the search region. While here this is for random parameters, this happens in real examples as well. This implements the parameters using the PETSc bag and thus they are all changeable on the commandline. Need to get some realistic parameters in place.

