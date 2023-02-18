// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro;

/** \brief A class for looping over all particles in a container_periodic or
 * container_periodic_poly class.
 *
 * Since the container_periodic and container_periodic_poly classes have a
 * fundamentally different memory organization, the regular loop classes cannot
 * be used with them. */
public class CLoopAllPeriodic extends CLoopBase {

    /** The constructor copies several necessary constants from the
     * base periodic container class.
     * \param[in] con the periodic container class to use. */
    public CLoopAllPeriodic(ContainerPeriodicBase con) {
        super(con);
        ey = con.ey;
        ez = con.ez;
        wy = con.wy;
        wz = con.wz;
        ijk0 = nx*(ey+con.oy*ez);
        inc2 = 2*nx*con.ey+1;
    }
    /** Sets the class to consider the first particle.
     * \return True if there is any particle to consider, false
     * otherwise. */
    public boolean start() {
        i=0;
        j=ey;
        k=ez;
        ijk=ijk0;
        q=0;
        while(co[ijk]==0) if(!next_block()) return false;
        return true;
    }
    /** Finds the next particle to test.
     * \return True if there is another particle, false if no more
     * particles are available. */
    public boolean inc() {
        q++;
        if(q>=co[ijk]) {
            q=0;
            do {
                if(!next_block()) return false;
            } while(co[ijk]==0);
        }
        return true;
    }

    /** The lower y index (inclusive) of the primary domain within
     * the block structure. */
    private final int ey;
    /** The lower y index (inclusive) of the primary domain within
     * the block structure. */
    private final int ez;
    /** The upper y index (exclusive) of the primary domain within
     * the block structure. */
    private final int wy;
    /** The upper z index (exclusive) of the primary domain within
     * the block structure. */
    private final int wz;
    /** The index of the (0,0,0) block within the block structure.
     */
    private final int ijk0;
    /** A value to increase ijk by when the z index is increased.
     */
    private final int inc2;
    /** Updates the internal variables to find the next
     * computational block with any particles.
     * \return True if another block is found, false if there are
     * no more blocks. */
    private boolean next_block() {
        i++;
        if(i==nx) {
            i=0;j++;
            if(j==wy) {
                j=ey;k++;
                if(k==wz) return false;
                ijk+=inc2;
            } else ijk++;
        } else ijk++;
        return true;
    }

}
