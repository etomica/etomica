// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro;

import java.util.Arrays;

/** \brief A class for storing an arbitrary number of particles, prior to setting
 * up a container geometry.
 *
 * The pre_container_base class can dynamically import and store an arbitrary
 * number of particles. Once the particles have been read in, an appropriate
 * container class can be set up with the optimal grid size, and the particles
 * can be transferred.
 *
 * The pre_container_base class is not intended for direct use, but forms the
 * base of the pre_container and pre_container_poly classes, that add routines
 * depending on whether particle radii need to be tracked or not. */
public class PreContainerBase {

    /** The minimum x coordinate of the container. */
    public final double ax;
    /** The maximum x coordinate of the container. */
    public final double bx;
    /** The minimum y coordinate of the container. */
    public final double ay;
    /** The maximum y coordinate of the container. */
    public final double by;
    /** The minimum z coordinate of the container. */
    public final double az;
    /** The maximum z coordinate of the container. */
    public final double bz;
    /** A boolean value that determines if the x coordinate in
     * periodic or not. */
    public final boolean xperiodic;
    /** A boolean value that determines if the y coordinate in
     * periodic or not. */
    public final boolean yperiodic;
    /** A boolean value that determines if the z coordinate in
     * periodic or not. */
    public final boolean zperiodic;

    /** Makes a guess at the optimal grid of blocks to use, computing in
     * a way that
     * \param[out] (nx,ny,nz) the number of blocks to use. */
    public void guess_optimal(int[] nx,int[] ny,int[] nz) {
        double dx=bx-ax,dy=by-ay,dz=bz-az;
        double ilscale=Math.pow(total_particles()/(Config.optimal_particles*dx*dy*dz),1/3.0);
        nx[0]=(int)(dx*ilscale+1);
        ny[0]=(int)(dy*ilscale+1);
        nz[0]=(int)(dz*ilscale+1);
    }

    /** The class constructor sets up the geometry of container, initializing the
     * minimum and maximum coordinates in each direction. It allocates an initial
     * chunk into which to store particle information.
     * \param[in] (ax_,bx_) the minimum and maximum x coordinates.
     * \param[in] (ay_,by_) the minimum and maximum y coordinates.
     * \param[in] (az_,bz_) the minimum and maximum z coordinates.
     * \param[in] (xperiodic_,yperiodic_,zperiodic_ ) flags setting whether the
     *                                                container is periodic in each
     *                                                coordinate direction.
     * \param[in] ps_ the number of floating point entries to store for each
     *                particle. */
    public PreContainerBase(double ax_,double bx_,double ay_,double by_,double az_,double bz_,boolean xperiodic_,boolean yperiodic_,boolean zperiodic_,int ps_) {
        ax = ax_;
        bx = bx_;
        ay = ay_;
        by = by_;
        az = az_;
        bz = bz_;
        xperiodic = xperiodic_;
        yperiodic = yperiodic_;
        zperiodic = zperiodic_;
        ps = ps_;
        int index_sz = Config.init_chunk_size;
        pre_id = new int[index_sz][];
        end_id = 0; // pre_id
        pre_id[end_id] = new int[Config.pre_container_chunk_size];
        pre_p = new double[index_sz][];
        end_p = 0;
        ch_id = 0; // pre_id[endid]
        e_id = ch_id + Config.pre_container_chunk_size;
        pre_p[end_p] = new double[ps*Config.pre_container_chunk_size];
        ch_p = 0; // pre_p[end_p]
    }
    /** Calculates and returns the total number of particles stored
     * within the class.
     * \return The number of particles. */
    public int total_particles() {
        return end_id*Config.pre_container_chunk_size+ch_id;
    }

    /** The number of doubles associated with a single particle
     * (three for the standard container, four when radius
     * information is stored). */
    protected final int ps;

    /** Allocates a new chunk of memory for storing particles. */
    protected void new_chunk() {
        end_id++;
        end_p++;
        if(end_id==pre_id.length) extend_chunk_index();
        ch_id=0;
        pre_id[end_id]=new int[Config.pre_container_chunk_size];
        e_id=ch_id+Config.pre_container_chunk_size;
        pre_p[end_p] = new double[ps*Config.pre_container_chunk_size];
        ch_p=0;
    }

    /** Extends the index of chunks. */
    protected void extend_chunk_index() {
        int index_sz = pre_id.length <<1;
        if(index_sz>Config.max_chunk_size)
            Common.voro_fatal_error("Absolute memory limit on chunk index reached",Config.Voropp.MEMORY_ERROR);
        if (Config.VOROPP_VERBOSE >= 2) {
            System.err.printf("Pre-container chunk index scaled up to %d\n", index_sz);
        }
        pre_id = Arrays.copyOf(pre_id, index_sz);
        pre_p = Arrays.copyOf(pre_p, index_sz);
    }

    /** A pointer to the chunk index to store the integer particle
     * IDs. */
    protected int[][] pre_id;
    /** A pointer to the last allocated integer ID chunk. */
    protected int end_id;
    /** A pointer to the next available slot on the current
     * particle ID chunk. */
    protected int ch_id;
    /** A pointer to the end of the current integer chunk. */
    protected int e_id;
    /** A pointer to the chunk index to store the floating point
     * information associated with particles. */
    protected double[][] pre_p;
    /** A pointer to the last allocated chunk of floating point
     * information. */
    protected int end_p;
    /** A pointer to the next available slot on the current
     * floating point chunk. */
    protected int ch_p;

}
