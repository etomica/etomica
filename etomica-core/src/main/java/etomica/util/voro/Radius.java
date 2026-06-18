/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.util.voro;

public interface Radius {
    void r_init(int ijk, int s);

    void r_prime(double rv);

    boolean r_ctest(double crs, double mrs);

    double r_cutoff(double lrs);

    double r_max_add(double rs);

    double r_current_sub(double rs, int ijk, int q);

    double r_scale(double rs, int ijk, int q);

    boolean r_scale_check(double[] rs, double mrs, int ijk, int q);
}
