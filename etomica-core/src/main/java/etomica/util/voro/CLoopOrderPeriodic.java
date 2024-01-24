// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro;

/** \brief Class for looping over all of the particles specified in a
 * pre-assembled particle_order class, for use with container_periodic classes.
 *
 * The particle_order class can be used to create a specific order of particles
 * within the container. This class can then loop over these particles in this
 * order. The class is particularly useful in cases where the ordering of the
 * output must match the ordering of particles as they were inserted into the
 * container. */
public class CLoopOrderPeriodic extends CLoopBase {

    /** A reference to the ordering class to use. */
    ParticleOrder vo;
    /** A pointer to the current position in the ordering class. */
    int cp;
    /** A pointer to the end position in the ordering class. */
    int op;
    /** The constructor copies several necessary constants from the
     * base class, and sets up a reference to the ordering class to
     * use.
     * \param[in] con the container class to use.
     * \param[in] vo_ the ordering class to use. */
    public CLoopOrderPeriodic(ContainerPeriodicBase con, ParticleOrder vo_) {
        super(con);
        vo = vo_;
        nx = con.nx;
        oxy = con.nx * con.oy;
    }
    /** Sets the class to consider the first particle.
     * \return True if there is any particle to consider, false
     * otherwise. */
    public boolean start() {
        cp = 0; // vo.o
        op=vo.op;
        if(cp!=op) {
            ijk = vo.o[cp];
            cp++;
            decode();
            q = vo.o[cp];
            cp++;
            return true;
        } else return false;
    }
    /** Finds the next particle to test.
     * \return True if there is another particle, false if no more
     * particles are available. */
    public boolean inc() {
        if(cp==op) return false;
        ijk = vo.o[cp];
        cp++;
        decode();
        q = vo.o[cp];
        cp++;
        return true;
    }

    /** The number of computational blocks in the x direction. */
    private final int nx;
    /** The number of computational blocks in a z-slice. */
    private final int oxy;
    /** Takes the current block index and computes indices in the
     * x, y, and z directions. */
    private void decode() {
        k=ijk/oxy;
        int ijkt=ijk-oxy*k;
        j=ijkt/nx;
        i=ijkt-j*nx;
    }


}
