/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.util.voro;

public class Config {

    // These constants set the initial memory allocation for the Voronoi cell
/** The initial memory allocation for the number of vertices. */
public static final int init_vertices=256;
/** The initial memory allocation for the maximum vertex order. */
public static final int init_vertex_order=64;
/** The initial memory allocation for the number of regular vertices of order
 * 3. */
public static final int init_3_vertices=256;
/** The initial memory allocation for the number of vertices of higher order.
 */
public static final int init_n_vertices=8;
/** The initial size for the delete stack. */
public static final int init_delete_size=256;
/** The initial size for the auxiliary delete stack. */
public static final int init_delete2_size=256;
/** The initial size for the extra search stack. */
public static final int init_xsearch_size=32;
/** The initial size for the wall pointer array. */
public static final int init_wall_size=32;
/** The default initial size for the ordering class. */
public static final int init_ordering_size=4096;
/** The initial size of the pre_container chunk index. */
public static final int init_chunk_size=256;

// If the initial memory is too small, the program dynamically allocates more.
// However, if the limits below are reached, then the program bails out.
/** The maximum memory allocation for the number of vertices. */
public static final int max_vertices=16777216;
/** The maximum memory allocation for the maximum vertex order. */
public static final int max_vertex_order=2048;
/** The maximum memory allocation for the any particular order of vertex. */
public static final int max_n_vertices=16777216;
/** The maximum size for the delete stack. */
public static final int max_delete_size=16777216;
/** The maximum size for the auxiliary delete stack. */
public static final int max_delete2_size=16777216;
/** The maximum size for the extra search stack. */
public static final int max_xsearch_size=16777216;
/** The maximum amount of particle memory allocated for a single region. */
public static final int max_particle_memory=16777216;
/** The maximum size for the wall pointer array. */
public static final int max_wall_size=2048;
/** The maximum size for the ordering class. */
public static final int max_ordering_size=67108864;
/** The maximum size for the pre_container chunk index. */
public static final int max_chunk_size=65536;

/** The chunk size in the pre_container classes. */
public static final int pre_container_chunk_size=1024;

/** Voro++ can print a number of different status and debugging messages to
 * notify the user of special behavior, and this macro sets the amount which
 * are displayed. At level 0, no messages are printed. At level 1, messages
 * about unusual cases during cell construction are printed, such as when the
 * plane routine bails out due to floating point problems. At level 2, general
 * messages about memory expansion are printed. At level 3, technical details
 * about memory management are printed. */
public static final int VOROPP_VERBOSE = 2;

/** If a point is within this distance of a cutting plane, then the code
 * assumes that point exactly lies on the plane. */
public static final double tolerance=1e-10*2.22045e-16;

public static final double big_tolerance_fac=20.;

public static double default_length=1000.;

/** A large number that is used in the computation. */
public static final double large_number=Double.MAX_VALUE;

/** A radius to use as a placeholder when no other information is available. */
public static final double default_radius=0.5;

/** The maximum number of shells of periodic images to test over. */
public static final int max_unit_voro_shells=10;

/** A guess for the optimal number of particles per block, used to set up the
 * container grid. */
public static final double optimal_particles=5.6;

/** If this is set to 1, then the code reports any instances of particles being
 * put outside of the container geometry. */
public static final boolean VOROPP_REPORT_OUT_OF_BOUNDS = false;

public enum Voropp {
    /** Voro++ returns this status code if there is a file-related error, such as
     * not being able to open file. */
    FILE_ERROR,
    /** Voro++ returns this status code if there is a memory allocation error, if
     * one of the safe memory limits is exceeded. */
    MEMORY_ERROR,
    /** Voro++ returns this status code if there is any type of internal error, if
     * it detects that representation of the Voronoi cell is inconsistent. This
     * status code will generally indicate a bug, and the developer should be
     * contacted. */
    INTERNAL_ERROR,
    /** Voro++ returns this status code if it could not interpret the command line
     * arguments passed to the command line utility. */
    CMD_LINE_ERROR
};

}
