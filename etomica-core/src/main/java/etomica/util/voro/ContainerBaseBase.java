/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.util.voro;

public abstract class ContainerBaseBase extends VoroBase implements Radius {
    /**
     * The maximum length squared that could be encountered in the
     * Voronoi cell calculation.
     */
    public double max_len_sq;

    /**
     * A two dimensional array holding particle positions. For the
     * derived container_poly class, this also holds particle
     * radii.
     */
    public double[][] p;
    /**
     * This array holds the numerical IDs of each particle in each
     * computational box.
     */
    public int[][] id;
    /**
     * This array holds the number of particles within each
     * computational box of the container.
     */
    public int[] co;
    /**
     * The amount of memory in the array structure for each
     * particle. This is set to 3 when the basic class is
     * initialized, so that the array holds (x,y,z) positions. If
     * the container class is initialized as part of the derived
     * class container_poly, then this is set to 4, to also hold
     * the particle radii.
     */
    public final int ps;
    public Radius radius;

    public ContainerBaseBase(int nx_, int ny_, int nz_, double boxx_, double boxy_, double boxz_, int ps_) {
        super(nx_, ny_, nz_, boxx_, boxy_, boxz_);
        ps = ps_;
    }

    public double r_current_sub(double rs, int ijk, int q) {
        return radius.r_current_sub(rs, ijk, q);
    }

    public void r_init(int ijk, int s) {
        radius.r_init(ijk, s);
    }

    public boolean r_scale_check(double[] rs, double mrs, int ijk, int q) {
        return radius.r_scale_check(rs, mrs, ijk, q);
    }

    public void r_prime(double rv) {
        radius.r_prime(rv);
    }

    public boolean r_ctest(double crs, double mrs) {
        return radius.r_ctest(crs, mrs);
    }

    public double r_cutoff(double lrs) {
        return radius.r_cutoff(lrs);
    }

    public double r_max_add(double rs) {
        return radius.r_max_add(rs);
    }

    public double r_scale(double rs, int ijk, int q) {
        return radius.r_scale(rs, ijk, q);
    }


    public abstract boolean initialize_voronoicell(VoronoiCellBase[] c, int ijk, int q, int ci, int cj, int ck,
                                                   int[] i, int[] j, int[] k, double[] x, double[] y, double[] z, int[] disp);

    public abstract void frac_pos(double x, double y, double z, double ci, double cj, double ck,
                                  double[] fx, double[] fy, double[] fz);

    public abstract int region_index(int ci,int cj,int ck,int ei,int ej,int ek,double[] qx,double[] qy,double[] qz,int[] disp);

    public abstract void initialize_search(int ci,int cj,int ck,int ijk,int[] i,int[] j,int[] k,int[] disp);

}