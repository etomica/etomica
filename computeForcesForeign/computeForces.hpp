void computeForces(double *xs, double *ys, double *zs, double *fxs, double *fys, double *fzs, int n, double boxSize);

void updatePreForces(double *xs, double *ys, double *zs,
                                double *vxs, double *vys, double *vzs,
                                double *fxs, double *fys, double *fzs,
                                int n, double rm, double timeStep);

extern "C" void updatePostForces(double *vxs, double *vys, double *vzs,
                                 double *fxs, double *fys, double *fzs,
                                 int n, double rm, double timeStep);