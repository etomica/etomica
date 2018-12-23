#include <cstddef>
#include <cstdio>
#include <cstring>
#include <immintrin.h>
#include <iostream>
#include <iomanip>

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

void printVec(const __m256d &vec) {
    double buf[4];
    _mm256_storeu_pd(buf, vec);
    printf("(%f %f %f %f)\n", buf[0], buf[1], buf[2], buf[3]);
}

void printVecHex(const __m256d &vec) {
    double buf[4];
    _mm256_storeu_pd(buf, vec);
    std::cout << std::hex;
    for(int i = 0; i < 4; i++) {
        std::cout << *reinterpret_cast<uint64_t*>(&buf[i]) << " ";
    }
    std::cout << std::endl;
}

__m256d nearestImageVec(__m256d dists, double dimComponent) {
    double half = dimComponent / 2;
    __m256d dMut = dists;
    __m256d vecDim = _mm256_set1_pd(dimComponent);
    __m256d vecHalf = _mm256_set1_pd(half);
    __m256d mask = _mm256_cmp_pd(dMut, vecHalf, _CMP_GT_OQ);

    //todo: while loop
    dMut = _mm256_blendv_pd(dMut, _mm256_sub_pd(dMut, vecDim), mask);

    vecHalf = _mm256_set1_pd(-half);
    mask = _mm256_cmp_pd(dMut, vecHalf, _CMP_LT_OQ);

    dMut = _mm256_blendv_pd(dMut, _mm256_add_pd(dMut, vecDim), mask);
    return dMut;
}

__m256d duVec(__m256d r2s) {
    __m256d s2 = _mm256_div_pd(_mm256_set1_pd(1.0), r2s);
    __m256d s6 = _mm256_mul_pd(_mm256_mul_pd(s2, s2), s2);
    __m256d du = _mm256_mul_pd(
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
        int j = i + 1;
        int nbrCount = 0;
        int nbrSize = n - j;
        int vectorizedSize = nbrSize & (-4);

        __m256d ix = _mm256_broadcast_sd(xs + i);
        __m256d iy = _mm256_broadcast_sd(ys + i);
        __m256d iz = _mm256_broadcast_sd(zs + i);

        for(; nbrCount < vectorizedSize; j += 4, nbrCount += 4) {
            __m256d jx = _mm256_loadu_pd(xs + j);
            __m256d jy = _mm256_loadu_pd(ys + j);
            __m256d jz = _mm256_loadu_pd(zs + j);

            __m256d dx = _mm256_sub_pd(jx, ix);
            __m256d dy = _mm256_sub_pd(jy, iy);
            __m256d dz = _mm256_sub_pd(jz, iz);
            // if (i == 0 && j == 21)
            // {
            //     printVec(dz);
            // }

            dx = nearestImageVec(dx, boxSize);
            dy = nearestImageVec(dy, boxSize);
            dz = nearestImageVec(dz, boxSize);
            // if (i == 0 && j == 21)
            // {
            //     printVec(dz);
            // }

            __m256d x2 = _mm256_mul_pd(dx, dx);
            __m256d y2 = _mm256_mul_pd(dy, dy);
            __m256d z2 = _mm256_mul_pd(dz, dz);

            __m256d r2s = _mm256_add_pd(x2, _mm256_add_pd(y2, z2));
            // if (i == 0 && j == 21)
            // {
            //     printVec(r2s);
            // }
            __m256d cutoffMask = _mm256_cmp_pd(r2s, _mm256_set1_pd(9.0), _CMP_GE_OS);
            // if (i == 0 && j == 21)
            // {
            //     printf("mask: ");
            //     printVecHex(cutoffMask);
            //     std::cout << std::hex << std::setfill('0') << std::setw(4) << _mm256_movemask_pd(cutoffMask) << std::endl;
            //     // printf("int: %d", _mm256_movemask_pd(cutoffMask));
            // }

            if(_mm256_movemask_pd(cutoffMask) == 0xF) {
                continue;
            }

            __m256d dus = duVec(r2s);
            // if (i == 0 && j == 21)
            // {
            //     printVec(dus);
            // }
            __m256d tmp = _mm256_div_pd(dus, r2s);
            dx = _mm256_mul_pd(dx, tmp);
            dy = _mm256_mul_pd(dy, tmp);
            dz = _mm256_mul_pd(dz, tmp);

            // ...

            double iForces[4];
            _mm256_storeu_pd(iForces, dx);
            fxs[i] += iForces[0] + iForces[1] + iForces[2] + iForces[3];
            _mm256_storeu_pd(iForces, dy);
            fys[i] += iForces[0] + iForces[1] + iForces[2] + iForces[3];
            _mm256_storeu_pd(iForces, dz);
            fzs[i] += iForces[0] + iForces[1] + iForces[2] + iForces[3];


            _mm256_storeu_pd(fxs + j, _mm256_sub_pd(_mm256_loadu_pd(fxs + j), dx));
            _mm256_storeu_pd(fys + j, _mm256_sub_pd(_mm256_loadu_pd(fys + j), dy));
            _mm256_storeu_pd(fzs + j, _mm256_sub_pd(_mm256_loadu_pd(fzs + j), dz));
        }

        for (; nbrCount < nbrSize; j++, nbrCount++) {
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

int main(int argc, char const *argv[])
{
    __m256d mask = _mm256_cmp_pd(_mm256_set_pd(1, 9, 3, 4), _mm256_set1_pd(2), _CMP_GT_OQ);
    printVec(mask);
    std::cout << std::hex << std::setfill('0') << std::setw(4) << _mm256_movemask_pd(mask) << std::endl;
    return 0;
}
