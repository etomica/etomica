#include <cstddef>
#include <cstdio>
#include <cstring>
#include <immintrin.h>

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
    memset(fxs, 0, n * sizeof(double));
    memset(fys, 0, n * sizeof(double));
    memset(fzs, 0, n * sizeof(double));
    
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

__m256d nearestImageVec(__m256d dists, double dimComponent) {
    double half = dimComponent / 2;
    auto dMut = dists;
    auto vecDim = _mm256_set1_pd(dimComponent);
    auto vecHalf = _mm256_set1_pd(half);
    auto mask = _mm256_cmp_pd(dMut, vecHalf, _CMP_GT_OQ);

    //todo: while loop
    dMut = _mm256_blendv_pd(dMut, _mm256_sub_pd(dMut, vecDim), mask);

    vecHalf = _mm256_set1_pd(-half);
    mask = _mm256_cmp_pd(dMut, vecHalf, _CMP_LT_OQ);

    dMut = _mm256_blendv_pd(dMut, _mm256_add_pd(dMut, vecDim), mask);
    return dMut;
}

__m256d duVec(__m256d r2s) {
    auto s2 = _mm256_div_pd(_mm256_set1_pd(1.0), r2s);
    auto s6 = _mm256_mul_pd(_mm256_mul_pd(s2, s2), s2);
    auto du = _mm256_mul_pd(
        _mm256_mul_pd(_mm256_set1_pd(-48.0), s6),
        _mm256_sub_pd(s6, _mm256_set1_pd(0.5)));

    return _mm256_blendv_pd(
        du,
        _mm256_set1_pd(0.0),
        _mm256_cmp_pd(r2s, _mm256_set1_pd(9.0), _CMP_GE_OQ));
}

extern "C" void computeForcesSimd(double *xs, double *ys, double *zs, double *fxs, double *fys, double *fzs, int n, double boxSize) {
    memset(fxs, 0, n * sizeof(double));
    memset(fys, 0, n * sizeof(double));
    memset(fzs, 0, n * sizeof(double));

    for(int i = 0; i < n; i++) {
        int j = i + i;
        int nbrCount = 0;
        int nbrSize = n - j;
        int vectorizedSize = nbrSize & (-4);

        __m256d ix = _mm256_loadu_pd(xs + i);
        __m256d iy = _mm256_loadu_pd(ys + i);
        __m256d iz = _mm256_loadu_pd(zs + i);

        for(; nbrCount < vectorizedSize; j += 4, nbrCount += 4) {
            __m256d dx = _mm256_sub_pd(_mm256_loadu_pd(xs + j), ix);
            __m256d dy = _mm256_sub_pd(_mm256_loadu_pd(ys + j), iy);
            __m256d dz = _mm256_sub_pd(_mm256_loadu_pd(zs + j), iz);

            dx = nearestImageVec(dx, boxSize);
            dy = nearestImageVec(dy, boxSize);
            dz = nearestImageVec(dz, boxSize);

            __m256d x2 = _mm256_mul_pd(dx, dx);
            __m256d y2 = _mm256_mul_pd(dy, dy);
            __m256d z2 = _mm256_mul_pd(dz, dz);

            __m256d r2s = _mm256_add_pd(x2, _mm256_add_pd(y2, z2));
            __m256d cutoffMask = _mm256_cmp_pd(r2s, _mm256_set1_pd(9.0), _CMP_GE_OQ);
            if(_mm256_movemask_pd(cutoffMask) == 0xff) {
                continue;
            }

            __m256d dus = duVec(r2s);
            __m256d tmp = _mm256_div_pd(dus, r2s);
            dx = _mm256_mul_pd(dx, tmp);
            dy = _mm256_mul_pd(dy, tmp);
            dz = _mm256_mul_pd(dz, tmp);

            // ...
        }
    }

}

extern "C" double computeEnergy(double *xs, double *ys, double *zs, int n, double boxSize) {
    double energy = 0;
    for(int i = 0; i < n; i++) {
        for(int j = i + 1; j < n; j++) {
            double dx = xs[j] - xs[i];
            double dy = ys[j] - ys[i];
            double dz = zs[j] - zs[i];

            dx = nearestImage(dx, boxSize);
            dy = nearestImage(dy, boxSize);
            dz = nearestImage(dz, boxSize);
            double r2 = dx * dx + dy * dy + dz * dz;

            if(r2 > 3 * 3) { continue; }

            double s2 = 1.0 / r2;
            double s6 = s2 * s2 * s2;
            double u = 4.0 * s6 * (s6 - 1.0);
            energy += u;
        }
    }
    return energy;
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