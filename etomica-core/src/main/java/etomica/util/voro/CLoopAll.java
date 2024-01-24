// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro;

/** \brief Class for looping over all of the particles in a container.
 *
 * This is one of the simplest loop classes, that scans the computational
 * blocks in order, and scans all the particles within each block in order. */
public class CLoopAll extends CLoopBase {

    /** The constructor copies several necessary constants from the
     * base container class.
     * \param[in] con the container class to use. */
    public CLoopAll(ContainerBase con) {
        super(con);
    }
    /** Sets the class to consider the first particle.
     * \return True if there is any particle to consider, false
     * otherwise. */
    public boolean start() {
        i=j=k=ijk=q=0;
        while(co[ijk]==0) {
            if(!next_block()) return false;
        }
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

    /** Updates the internal variables to find the next
     * computational block with any particles.
     * \return True if another block is found, false if there are
     * no more blocks. */
    private boolean next_block() {
        ijk++;
        i++;
        if(i==nx) {
            i=0;j++;
            if(j==ny) {
                j=0;k++;
                if(ijk==nxyz) return false;
            }
        }
        return true;
    }

}
