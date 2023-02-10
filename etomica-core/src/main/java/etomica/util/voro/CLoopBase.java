// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro;

/** \brief Base class for looping over particles in a container.
 *
 * This class forms the base of all classes that can loop over a subset of
 * particles in a contaner in some order. When initialized, it stores constants
 * about the corresponding container geometry. It also contains a number of
 * routines for interrogating which particle currently being considered by the
 * loop, which are common between all of the derived classes. */
public abstract class CLoopBase {

    /** The number of blocks in the x direction. */
    public final int nx;
    /** The number of blocks in the y direction. */
    public final int ny;
    /** The number of blocks in the z direction. */
    public final int nz;
    /** A constant, set to the value of nx multiplied by ny, which
     * is used in the routines that step through blocks in
     * sequence. */
    public final int nxy;
    /** A constant, set to the value of nx*ny*nz, which is used in
     * the routines that step through blocks in sequence. */
    public final int nxyz;
    /** The number of floating point numbers per particle in the
     * associated container data structure. */
    public final int ps;
    /** A pointer to the particle position information in the
     * associated container data structure. */
    public double [][] p;
    /** A pointer to the particle ID information in the associated
     * container data structure. */
    public int[][] id;
    /** A pointer to the particle counts in the associated
     * container data structure. */
    public int[] co;
    /** The current x-index of the block under consideration by the
     * loop. */
    public int i;
    /** The current y-index of the block under consideration by the
     * loop. */
    public int j;
    /** The current z-index of the block under consideration by the
     * loop. */
    public int k;
    /** The current index of the block under consideration by the
     * loop. */
    public int ijk;
    /** The index of the particle under consideration within the current
     * block. */
    public int q;
    /** The constructor copies several necessary constants from the
     * base container class.
     * \param[in] con the container class to use. */
    public CLoopBase(ContainerBaseBase con) {
        nx = con.nx;
        ny = con.ny;
        nz = con.nz;
        nxy = con.nxy;
        nxyz = con.nxyz;
        ps = con.ps;
        p = con.p;
        id = con.id;
        co = con.co;
    }
    /** Returns the position vector of the particle currently being
     * considered by the loop.
     * \param[out] (x,y,z) the position vector of the particle. */
    public void pos(double[] x,double[] y,double[] z) {
        int pp = ps*q; //p[ijk]
        x[0]=p[ijk][pp];
        y[0]=p[ijk][pp+1];
        z[0]=p[ijk][pp+2];
    }
    /** Returns the ID, position vector, and radius of the particle
     * currently being considered by the loop.
     * \param[out] pid the particle ID.
     * \param[out] (x,y,z) the position vector of the particle.
     * \param[out] r the radius of the particle. If no radius
     * 		 information is available the default radius
     * 		 value is returned. */
    public void pos(int[] pid,double[] x,double[] y,double[] z,double[] r) {
        pid[0]=id[ijk][q];
        pos(x,y,z);
        r[0] = ps==3?Config.default_radius:p[ijk][ps*q+3];
    }
    /** Returns the x position of the particle currently being
     * considered by the loop. */
    public double x() {return p[ijk][ps*q];}
    /** Returns the y position of the particle currently being
     * considered by the loop. */
    public double y() {return p[ijk][ps*q+1];}
    /** Returns the z position of the particle currently being
     * considered by the loop. */
    public double z() {return p[ijk][ps*q+2];}
    /** Returns the ID of the particle currently being considered
     * by the loop. */
    public int pid() {return id[ijk][q];}

    public abstract boolean start();
    public abstract boolean inc();

}
