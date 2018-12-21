#include <cstddef>
#include <cstdio>
#include <cstring>

double du(double r2) {
    double s2 = 1.0 / r2;
    double s6 = s2 * s2 * s2;
    return -48.0 * s6 * (s6 - 0.5);
}

double nearestImage(double dr, double boxSize) {
    double half = boxSize / 2;
    double image = dr;
    while(image > half) {
        image -= boxSize;
    }

    while(image < -half) {
        image += boxSize;
    }

    return image;
}

extern "C" void computeForces(double *xs, double *ys, double *zs, double *fxs, double *fys, double *fzs, int n, double boxSize) {
    memset(fxs, 0, n);
    memset(fys, 0, n);
    memset(fzs, 0, n);
    
    for(size_t i = 0; i < n; i++) {
        
        for(size_t j = i + 1; j < n; j++) {
            double dx = xs[j] - xs[i];
            double dy = ys[j] - ys[i];
            double dz = zs[j] - zs[i];

            dx = nearestImage(dx, boxSize);
            dy = nearestImage(dy, boxSize);
            dz = nearestImage(dz, boxSize);

            double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 > 3 * 3) { continue; }

            double div = du(r2) / r2;
            dx *= div;
            dy *= div;
            dz *= div;
            
            fxs[i] += dx;
            fys[i] += dy;
            fzs[i] += dz;
            
            fxs[j] -= dx;
            fys[j] -= dy;
            fzs[j] -= dz;
        }
        // if(i == 0) {
        //     printf("%f %f %f\n", xs[i], ys[i], zs[i]);
        //     printf("force: %f %f %f\n", fxs[i], fys[i], fzs[i]);
        // }
        
    }
}

extern "C" void updatePreForces(double *xs, double *ys, double *zs,
                                double *vxs, double *vys, double *vzs,
                                double *fxs, double *fys, double *fzs,
                                int n, double rm, double timeStep)
{
    for(size_t i = 0; i < n; i++) {
        vxs[i] += 0.5 * rm * timeStep * fxs[i];
        vys[i] += 0.5 * rm * timeStep * fys[i];
        vzs[i] += 0.5 * rm * timeStep * fzs[i];

        xs[i] += timeStep * vxs[i];
        ys[i] += timeStep * vys[i];
        zs[i] += timeStep * vzs[i];
    }
}

extern "C" void updatePostForces(double *vxs, double *vys, double *vzs,
                                 double *fxs, double *fys, double *fzs,
                                 int n, double rm, double timeStep)
{
    for(size_t i = 0; i < n; i++) {
        vxs[i] += 0.5 * rm * timeStep * fxs[i];
        vys[i] += 0.5 * rm * timeStep * fys[i];
        vzs[i] += 0.5 * rm * timeStep * fzs[i];
    }
}