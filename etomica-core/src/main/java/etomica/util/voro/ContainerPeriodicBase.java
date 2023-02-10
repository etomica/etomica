// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro;

import etomica.util.collections.DoubleArrayList;
import etomica.util.collections.IntArrayList;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Arrays;

/** \brief Class for representing a particle system in a 3D periodic
 * non-orthogonal periodic domain.
 *
 * This class represents a particle system in a three-dimensional
 * non-orthogonal periodic domain. The domain is defined by three periodicity
 * vectors (bx,0,0), (bxy,by,0), and (bxz,byz,bz) that represent a
 * parallelepiped. Internally, the class stores particles in the box 0<x<bx,
 * 0<y<by, 0<z<bz, and constructs periodic images of particles that displaced
 * by the three periodicity vectors when they are necessary for the
 * computation. The internal memory structure for this class is significantly
 * different from the container_base class in order to handle the dynamic
 * construction of these periodic images.
 *
 * The class is derived from the unitcell class, which encapsulates information
 * about the domain geometry, and the voro_base class, which encapsulates
 * information about the underlying computational grid. */
public class ContainerPeriodicBase extends ContainerBaseBase {

    // unitcell bits
    /** The x coordinate of the first vector defining the periodic
	 * domain. */
	public final double bx;
    /** The x coordinate of the second vector defining the periodic
	 * domain. */
    public final double bxy;
    /** The y coordinate of the second vector defining the periodic
	 * domain. */
	public final double by;
	/** The x coordinate of the third vector defining the periodic
	 * domain. */
	public final double bxz;
	/** The y coordinate of the third vector defining the periodic
	 * domain. */
	public final double byz;
	/** The z coordinate of the third vector defining the periodic
	 * domain. */
	public final double bz;
	/** The computed unit Voronoi cell corresponding the given
	 * 3D non-rectangular periodic domain geometry. */
	public VoronoiCell unit_voro;
	/** Draws an outline of the domain in Gnuplot format.
	 * \param[in] filename the filename to write to. */
	public  void draw_domain_gnuplot(String filename) {
        try {
            FileOutputStream fp = new FileOutputStream(filename);
            draw_domain_gnuplot(fp);
            fp.close();
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
	}
    public void draw_domain_gnuplot() {
        draw_domain_gnuplot(System.out);
    }

    /** Draws the periodic domain in gnuplot format.
     * \param[in] fp the file handle to write to. */
	void draw_domain_gnuplot(OutputStream fp) {
        try {
            fp.write(String.format("0 0 0\n%g 0 0\n%g %g 0\n%g %g 0\n", bx, bx + bxy, by, bxy, by).getBytes());
            fp.write(String.format("%g %g %g\n%g %g %g\n%g %g %g\n%g %g %g\n", bxy + bxz, by + byz, bz, bx + bxy + bxz, by + byz, bz, bx + bxz, byz, bz, bxz, byz, bz).getBytes());
            fp.write(String.format("0 0 0\n%g %g 0\n\n%g %g %g\n%g %g %g\n\n", bxy, by, bxz, byz, bz, bxy + bxz, by + byz, bz).getBytes());
            fp.write(String.format("%g 0 0\n%g %g %g\n\n%g %g 0\n%g %g %g\n\n", bx, bx + bxz, byz, bz, bx + bxy, by, bx + bxy + bxz, by + byz, bz).getBytes());
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }
	/** Draws an outline of the domain in Gnuplot format.
	 * \param[in] filename the filename to write to. */
	void draw_domain_pov(String filename) {
        try {
            FileOutputStream fp = new FileOutputStream(filename);
            draw_domain_pov(fp);
            fp.close();
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
	}
    public void draw_domain_pov() {
        draw_domain_pov(System.out);
    }

    /** Draws the periodic domain in POV-Ray format.
     * \param[in] fp the file handle to write to. */
	void draw_domain_pov(OutputStream fp) {
        try {
            fp.write(String.format("cylinder{0,0,0>,<%g,0,0>,rr}\n" +
                    "cylinder{<%g,%g,0>,<%g,%g,0>,rr}\n", bx, bxy, by, bx + bxy, by).getBytes());
            fp.write(String.format("cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n" +
                    "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n", bxz, byz, bz, bx + bxz, byz, bz, bxy + bxz, by + byz, bz, bx + bxy + bxz, by + byz, bz).getBytes());
            fp.write(String.format("cylinder{<0,0,0>,<%g,%g,0>,rr}\n" +
                    "cylinder{<%g,0,0>,<%g,%g,0>,rr}\n", bxy, by, bx, bx + bxy, by).getBytes());
            fp.write(String.format("cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n" +
                    "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n", bxz, byz, bz, bxy + bxz, by + byz, bz, bx + bxz, byz, bz, bx + bxy + bxz, by + byz, bz).getBytes());
            fp.write(String.format("cylinder{<0,0,0>,<%g,%g,%g>,rr}\n" +
                    "cylinder{<%g,0,0>,<%g,%g,%g>,rr}\n", bxz, byz, bz, bx, bx + bxz, byz, bz).getBytes());
            fp.write(String.format("cylinder{<%g,%g,0>,<%g,%g,%g>,rr}\n" +
                    "cylinder{<%g,%g,0>,<%g,%g,%g>,rr}\n", bxy, by, bxy + bxz, by + byz, bz, bx + bxy, by, bx + by + bxz, by + byz, bz).getBytes());
            fp.write(String.format("sphere{<0,0,0>,rr}\nsphere{<%g,0,0>,rr}\n" +
                    "sphere{<%g,%g,0>,rr}\nsphere{<%g,%g,0>,rr}\n", bx, bxy, by, bx + bxy, by).getBytes());
            fp.write(String.format("sphere{<%g,%g,%g>,rr}\nsphere{<%g,%g,%g>,rr}\n" +
                    "sphere{<%g,%g,%g>,rr}\nsphere{<%g,%g,%g>,rr}\n", bxz, byz, bz, bx + bxz, byz, bz, bxy + bxz, by + byz, bz, bx + bxy + bxz, by + byz, bz).getBytes());
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }

    /** Calculates whether the unit Voronoi cell intersects a given periodic image
     * of the domain.
     * \param[in] (dx,dy,dz) the displacement of the periodic image.
     * \param[out] vol the proportion of the unit cell volume within this image,
     *                 only computed in the case that the two intersect.
     * \return True if they intersect, false otherwise. */
	public boolean intersects_image(double dx,double dy,double dz,double[] vol) {
    	final double bxinv=1/bx,byinv=1/by,bzinv=1/bz,ivol=bxinv*byinv*bzinv;
        VoronoiCell c = new VoronoiCell(unit_voro);
        dx*=2;dy*=2;dz*=2;
        if(!c.plane(0,0,bzinv,dz+1)) return false;
        if(!c.plane(0,0,-bzinv,-dz+1)) return false;
        if(!c.plane(0,byinv,-byz*byinv*bzinv,dy+1)) return false;
        if(!c.plane(0,-byinv,byz*byinv*bzinv,-dy+1)) return false;
        if(!c.plane(bxinv,-bxy*bxinv*byinv,(bxy*byz-by*bxz)*ivol,dx+1)) return false;
        if(!c.plane(-bxinv,bxy*bxinv*byinv,(-bxy*byz+by*bxz)*ivol,-dx+1)) return false;
        vol[0]=c.volume()*ivol;
        return true;
    }
    public static class IntQueue {
        public int[] storage;
        public int front, end, size;
        public IntQueue() {
            storage = new int[10];
            size = front = end = 0;
        }

        public void push(int x) {
            if (size == storage.length) {
                int[] newStorage = new int[size*2];
                System.arraycopy(storage, front, newStorage, 0, size-front);
                if (front>0) System.arraycopy(storage, 0, newStorage, size-front, front);
                front = 0;
                end = size-1;
                storage = newStorage;
            }
            end = (end+1)%size;
            storage[end] = x;
            size++;
        }

        public int front() {
            if (size == 0) throw new RuntimeException("empty");
            return storage[front];
        }

        public void pop() {
            if (size == 0) throw new RuntimeException("empty");
            front = (front+1)%size;
            size--;
        }

        public boolean empty() {
            return size>0;
        }
    }

    /** Computes a list of periodic domain images that intersect the unit Voronoi cell.
     * \param[out] vi a vector containing triplets (i,j,k) corresponding to domain
     *                images that intersect the unit Voronoi cell, when it is
     *                centered in the middle of the primary domain.
     * \param[out] vd a vector containing the fraction of the Voronoi cell volume
     *                within each corresponding image listed in vi. */
	public void images(IntArrayList vi, DoubleArrayList vd) {
	    int ms2=Config.max_unit_voro_shells*2+1;
        int mss=ms2*ms2*ms2;
        boolean[] a=new boolean[mss];
        int ac=Config.max_unit_voro_shells*(1+ms2*(1+ms2));  // a
        int i,j,k;
        double[] vol = new double[1];

        // Initialize mask
        Arrays.fill(a, 0, ac, true);
        Arrays.fill(a, ac+1, mss, true);

        // Set up the queue and add (0,0,0) image to it
        IntQueue q = new IntQueue();
        q.push(0);q.push(0);q.push(0);

        while(!q.empty()) {

            // Read the next entry on the queue
            i=q.front();q.pop();
            j=q.front();q.pop();
            k=q.front();q.pop();

            // Check intersection of this image
            if(intersects_image(i,j,k,vol)) {

                // Add this entry to the output vectors
                vi.add(i);
                vi.add(j);
                vi.add(k);
                vd.add(vol[0]);

                // Add neighbors to the queue if they have not been
                // tested
                int ap=ac+i+ms2*(j+ms2*k);
                if(k>-Config.max_unit_voro_shells&&a[ap-ms2*ms2]) {q.push(i);q.push(j);q.push(k-1);a[ap-ms2*ms2]=false;}
                if(j>-Config.max_unit_voro_shells&&a[ap-ms2]) {q.push(i);q.push(j-1);q.push(k);a[ap-ms2]=false;}
                if(i>-Config.max_unit_voro_shells&&a[ap-1]) {q.push(i-1);q.push(j);q.push(k);a[ap-1]=false;}
                if(i<Config.max_unit_voro_shells&&a[ap+1]) {q.push(i+1);q.push(j);q.push(k);a[ap+1]=false;}
                if(j<Config.max_unit_voro_shells&&a[ap+ms2]) {q.push(i);q.push(j+1);q.push(k);a[ap+ms2]=false;}
                if(k<Config.max_unit_voro_shells&&a[ap+ms2*ms2]) {q.push(i);q.push(j);q.push(k+1);a[ap+ms2*ms2]=false;}
            }
        }
    }

	/** The maximum y-coordinate that could possibly cut the
	 * computed unit Voronoi cell. */
	protected double max_uv_y;
	/** The maximum z-coordinate that could possibly cut the
	 * computed unit Voronoi cell. */
	protected double max_uv_z;

    /** Applies a pair of opposing plane cuts from a periodic image point
     * to the unit Voronoi cell.
     * \param[in] (i,j,k) the index of the periodic image to consider. */
    private void unit_voro_apply(int i,int j,int k) {
        double x=i*bx+j*bxy+k*bxz,y=j*by+k*byz,z=k*bz;
        unit_voro.plane(x,y,z);
        unit_voro.plane(-x,-y,-z);
    }

    /** Tests to see if a shell of periodic images could possibly cut the periodic
     * unit cell.
     * \param[in] l the index of the shell to consider.
     * \return True if a point in the shell cuts the cell, false otherwise. */
    private boolean unit_voro_intersect(int l) {
        int i,j;
        if(unit_voro_test(l,0,0)) return true;
        for(i=1;i<l;i++) {
            if(unit_voro_test(l,i,0)) return true;
            if(unit_voro_test(-l,i,0)) return true;
        }
        for(i=-l;i<=l;i++) if(unit_voro_test(i,l,0)) return true;
        for(i=1;i<l;i++) for(j=-l+1;j<=l;j++) {
            if(unit_voro_test(l,j,i)) return true;
            if(unit_voro_test(-j,l,i)) return true;
            if(unit_voro_test(-l,-j,i)) return true;
            if(unit_voro_test(j,-l,i)) return true;
        }
        for(i=-l;i<=l;i++) for(j=-l;j<=l;j++) if(unit_voro_test(i,j,l)) return true;
        return false;
    }

    /** Tests to see if a plane cut from a particular periodic image will cut th
     * unit Voronoi cell.
     * \param[in] (i,j,k) the index of the periodic image to consider.
     * \return True if the image cuts the cell, false otherwise. */
    private boolean unit_voro_test(int i,int j,int k) {
        double x=i*bx+j*bxy+k*bxz,y=j*by+k*byz,z=k*bz;
        double rsq=x*x+y*y+z*z;
        return unit_voro.plane_intersects(x,y,z,rsq);
    }

    // done with unit cell methods

    /** The lower y index (inclusive) of the primary domain within
     * the block structure. */
    public int ey;
    /** The lower z index (inclusive) of the primary domain within
     * the block structure. */
    public int ez;
    /** The upper y index (exclusive) of the primary domain within
     * the block structure. */
    public int wy;
    /** The upper z index (exclusive) of the primary domain within
     * the block structure. */
    public int wz;
    /** The total size of the block structure (including images) in
     * the y direction. */
    public int oy;
    /** The total size of the block structure (including images) in
     * the z direction. */
    public int oz;
    /** The total number of blocks. */
    public int oxyz;
    /** This array holds the maximum amount of particle memory for
     * each computational box of the container. If the number of
     * particles in a particular box ever approaches this limit,
     * more is allocated using the add_particle_memory() function.
     */
    public int[] mem;
    /** An array holding information about periodic image
     * construction at a given location. */
    public byte[] img;
    /** The initial amount of memory to allocate for particles
     * for each block. */
    public final int init_mem;

    /** Initializes the unit cell class for a particular non-orthogonal periodic
     * geometry, corresponding to a parallelepiped with sides given by three
     * vectors. The class constructs the unit Voronoi cell corresponding to this
     * geometry.
     * \param[in] (bx_) The x coordinate of the first unit vector.
     * \param[in] (bxy_,by_) The x and y coordinates of the second unit vector.
     * \param[in] (bxz_,byz_,bz_) The x, y, and z coordinates of the third unit
     *                            vector. */

    /** The class constructor sets up the geometry of container, initializing the
     * minimum and maximum coordinates in each direction, and setting whether each
     * direction is periodic or not. It divides the container into a rectangular
     * grid of blocks, and allocates memory for each of these for storing particle
     * positions and IDs.
     * \param[in] (bx_) The x coordinate of the first unit vector.
     * \param[in] (bxy_,by_) The x and y coordinates of the second unit vector.
     * \param[in] (bxz_,byz_,bz_) The x, y, and z coordinates of the third unit
     *                            vector.
     * \param[in] (nx_,ny_,nz_) the number of grid blocks in each of the three
     *                       coordinate directions.
     * \param[in] init_mem_ the initial memory allocation for each block.
     * \param[in] ps_ the number of floating point entries to store for each
     *                particle. */
    public ContainerPeriodicBase(double bx_,double bxy_,double by_,double bxz_,double byz_,double bz_,
                            int nx_,int ny_,int nz_,int init_mem_,int ps_) {
        super(nx_, ny_, nz_, bx_, by_, bz_, ps_);
        bx = bx_;
        bxy = bxy_;
        by = by_;
        bxz = bxz_;
        byz = byz_;
        bz = bz_;
        unit_voro = new VoronoiCell(Config.max_unit_voro_shells*Config.max_unit_voro_shells*4*(bx*bx+by*by+bz*bz));

        int l=1;

        // Initialize the Voronoi cell to be a very large rectangular box
	    final double ucx=Config.max_unit_voro_shells*bx;
        final double ucy=Config.max_unit_voro_shells*by;
        final double ucz=Config.max_unit_voro_shells*bz;
        unit_voro.init(-ucx,ucx,-ucy,ucy,-ucz,ucz);

        boolean success = false;
        // Repeatedly cut the cell by shells of periodic image particles
        while(l<2*Config.max_unit_voro_shells) {

            // Check to see if any of the planes from the current shell
            // will cut the cell
            if(unit_voro_intersect(l)) {

                // If they do, apply the plane cuts from the current
                // shell
                unit_voro_apply(l,0,0);
                for(int i=1;i<l;i++) {
                    unit_voro_apply(l,i,0);
                    unit_voro_apply(-l,i,0);
                }
                for(int i=-l;i<=l;i++) unit_voro_apply(i,l,0);
                for(int i=1;i<l;i++) for(int j=-l+1;j<=l;j++) {
                    unit_voro_apply(l,j,i);
                    unit_voro_apply(-j,l,i);
                    unit_voro_apply(-l,-j,i);
                    unit_voro_apply(j,-l,i);
                }
                for(int i=-l;i<=l;i++) for(int j=-l;j<=l;j++) unit_voro_apply(i,j,l);
            } else {

                // Calculate a bound on the maximum y and z coordinates
                // that could possibly cut the cell. This is based upon
                // a geometric result that particles with z>l can't cut
                // a cell lying within the paraboloid
                // z<=(l*l-x*x-y*y)/(2*l). It is always a tighter bound
                // than the one based on computing the maximum radius
                // of a Voronoi cell vertex.
                max_uv_y=max_uv_z=0;
                double[] pts=unit_voro.pts;
                for (int pp=0; pp<4*unit_voro.p; pp+=4) {
                    double x=pts[pp];
                    double y=pts[pp+1];
                    double z = pts[pp+2];
                    double q=Math.sqrt(x*x+y*y+z*z);
                    if(y+q>max_uv_y) max_uv_y=y+q;
                    if(z+q>max_uv_z) max_uv_z=z+q;
                }
                max_uv_z*=0.5;
                max_uv_y*=0.5;
                success = true;
                break;
            }
            l++;
        }

        // If the routine makes it here, then the unit cell still hasn't been
        // completely bounded by the plane cuts. Give the memory error code,
        // because this is mainly a case of hitting a safe limit, than any
        // inherent problem.
        if (!success) Common.voro_fatal_error("Periodic cell computation failed",Config.Voropp.MEMORY_ERROR);

        // done with unitcell
        max_len_sq = unit_voro.max_radius_squared();
        ey = (int)(max_uv_y*ysp+1);
        ez = (int)(max_uv_z*zsp+1);
        wy = ny+ey;
        wz = nz+ez;
        oy = ny+2*ey;
        oz = nz+2*ez;
        oxyz = nx*oy*oz;
        id = new int[oxyz][];
        p = new double[oxyz][];
        co = new int[oxyz];
        mem = new int[oxyz];
        img = new byte[oxyz];
        init_mem = init_mem_;

        // Set up memory for the blocks in the primary domain
        for(int k=ez;k<wz;k++) for(int j=ey;j<wy;j++) for(int i=0;i<nx;i++) {
            l=i+nx*(j+oy*k);
            mem[l]=init_mem;
            id[l]=new int[init_mem];
            p[l]=new double[ps*init_mem];
        }
    }

    /** Prints all particles in the container, including those that
     * have been constructed in image blocks. */
    public void print_all_particles() {
        int ijk,q;
        for(ijk=0;ijk<oxyz;ijk++) for(q=0;q<co[ijk];q++)
            System.out.printf("%d %g %g %g\n",id[ijk][q],p[ijk][ps*q],p[ijk][ps*q+1],p[ijk][ps*q+2]);
    }

    /** Outputs the a list of all the container regions along with the number of
     * particles stored within each. */
    void region_count() {
        for(int k=0,cop=0;k<nz;k++) for(int j=0;j<ny;j++) for(int i=0;i<nx;i++)
            System.out.printf("Region (%d,%d,%d): %d particles\n",i,j,k,co[cop]);
    }
    /** Initializes the Voronoi cell prior to a compute_cell
     * operation for a specific particle being carried out by a
     * voro_compute class. The cell is initialized to be the
     * pre-computed unit Voronoi cell based on planes formed by
     * periodic images of the particle.
     * \param[in,out] c a reference to a voronoicell object.
     * \param[in] ijk the block that the particle is within.
     * \param[in] q the index of the particle within its block.
     * \param[in] (ci,cj,ck) the coordinates of the block in the
     * 			 container coordinate system.
     * \param[out] (i,j,k) the coordinates of the test block
     * 		       relative to the voro_compute
     * 		       coordinate system.
     * \param[out] (x,y,z) the position of the particle.
     * \param[out] disp a block displacement used internally by the
     *		    compute_cell routine.
     * \return False if the plane cuts applied by walls completely
     * removed the cell, true otherwise. */
    public boolean initialize_voronoicell(VoronoiCellBase[] c,int ijk,int q,int ci,int cj,int ck,int[] i,int[] j,int[] k,double[] x,double[] y,double[] z,int[] disp) {
        c[0]=unit_voro;
        int pp = ps*q;  // p[ijk]
        x[0] = p[ijk][pp];
        y[0] = p[ijk][pp+1];
        z[0] = p[ijk][pp+2];
        i[0]=nx;
        j[0]=ey;
        k[0]=ez;
        return true;
    }
    /** Initializes parameters for a find_voronoi_cell call within
     * the voro_compute template.
     * \param[in] (ci,cj,ck) the coordinates of the test block in
     * 			 the container coordinate system.
     * \param[in] ijk the index of the test block
     * \param[out] (i,j,k) the coordinates of the test block
     * 		       relative to the voro_compute
     * 		       coordinate system.
     * \param[out] disp a block displacement used internally by the
     *		    find_voronoi_cell routine (but not needed
     *		    in this instance.) */
    public void initialize_search(int ci,int cj,int ck,int ijk,int[] i,int[] j,int[] k,int[] disp) {
        i[0]=nx;
        j[0]=ey;
        k[0]=ez;
    }
    /** Returns the position of a particle currently being computed
     * relative to the computational block that it is within. It is
     * used to select the optimal worklist entry to use.
     * \param[in] (x,y,z) the position of the particle.
     * \param[in] (ci,cj,ck) the block that the particle is within.
     * \param[out] (fx,fy,fz) the position relative to the block.
     */
    public void frac_pos(double x,double y,double z,double ci,double cj,double ck,double[] fx,double[] fy,double[] fz) {
        fx[0]=x-boxx*ci;
        fy[0]=y-boxy*(cj-ey);
        fz[0]=z-boxz*(ck-ez);
    }
    /** Calculates the index of block in the container structure
     * corresponding to given coordinates.
     * \param[in] (ci,cj,ck) the coordinates of the original block
     * 			 in the current computation, relative
     * 			 to the container coordinate system.
     * \param[in] (ei,ej,ek) the displacement of the current block
     * 			 from the original block.
     * \param[in,out] (qx,qy,qz) the periodic displacement that
     * 			     must be added to the particles
     * 			     within the computed block.
     * \param[in] disp a block displacement used internally by the
     * 		    find_voronoi_cell and compute_cell routines
     * 		    (but not needed in this instance.)
     * \return The block index. */
    public int region_index(int ci,int cj,int ck,int ei,int ej,int ek,double[] qx,double[] qy,double[] qz,int[] disp) {
        int qi=ci+(ei-nx),qj=cj+(ej-ey),qk=ck+(ek-ez);
        int iv = step_div(qi,nx);
        if(iv!=0) {
            qx[0]=iv*bx;
            qi-=nx*iv;
        }
        else {
            qx[0]=0;
        }
        create_periodic_image(qi,qj,qk);
        return qi+nx*(qj+oy*qk);
    }

    /** This routine creates all periodic images of the particles. It is meant for
     * diagnostic purposes only, since usually periodic images are dynamically
     * created in when they are referenced. */
    public void create_all_images() {
        for(int k=0;k<oz;k++) for(int j=0;j<oy;j++) for(int i=0;i<nx;i++) create_periodic_image(i,j,k);
    }

    /** Checks that the particles within each block lie within that block's bounds.
     * This is useful for diagnosing problems with periodic image computation. */
    public void check_compartmentalized() {
        int l;
        for(int k=l=0;k<oz;k++) for(int j=0;j<oy;j++) for(int i=0;i<nx;i++,l++) if(mem[l]>0) {

            // Compute the block's bounds, adding in a small tolerance
            double mix=i*boxx-Config.tolerance;
            double max=mix+boxx+Config.tolerance;
            double miy=(j-ey)*boxy-Config.tolerance;
            double may=miy+boxy+Config.tolerance;
            double miz=(k-ez)*boxz-Config.tolerance;
            double maz=miz+boxz+Config.tolerance;

            // Print entries for any particles that lie outside the block's
            // bounds
            double[] pp = p[l];
            for(int c=0;c<co[l];c++) if(pp[c*ps]<mix||pp[c*ps]>max||pp[c*ps+1]<miy||pp[c*ps+1]>may||pp[c*ps+2]<miz||pp[c*ps+2]>maz)
            System.out.printf("%d %d %d %d %f %f %f %f %f %f %f %f %f\n",
                    id[l][c],i,j,k,pp[c*ps],pp[c*ps+1],pp[c*ps+2],mix,max,miy,may,miz,maz);
        }
    }

    /** Increase memory for a particular region.
     * \param[in] i the index of the region to reallocate. */
    protected void add_particle_memory(int i) {

        // Handle the case when no memory has been allocated for this block
        if(mem[i]==0) {
            mem[i]=init_mem;
            id[i]=new int[init_mem];
            p[i]=new double[ps*init_mem];
            return;
        }

        // Otherwise, double the memory allocation for this block. Carry out a
        // check on the memory allocation size, and print a status message if
        // requested.
        int l;
        int nmem = mem[i]<<1;
        if(nmem>Config.max_particle_memory)
            Common.voro_fatal_error("Absolute maximum memory allocation exceeded",Config.Voropp.MEMORY_ERROR);
        if (Config.VOROPP_VERBOSE >= 3) {
            System.err.printf("Particle memory in region %d scaled up to %d\n", i, nmem);
        }

        // Allocate new memory and copy in the contents of the old arrays
        id[i] = Arrays.copyOf(id[i], nmem);
        p[i] = Arrays.copyOf(p[i], ps*nmem);

        // Update pointers and delete old arrays
        mem[i]=nmem;
    }

    /** Takes a particle position vector and computes the region index into which
     * it should be stored. If the container is periodic, then the routine also
     * maps the particle position to ensure it is in the primary domain. If the
     * container is not periodic, the routine bails out.
     * \param[out] ijk the region index.
     * \param[in,out] (x,y,z) the particle position, remapped into the primary
     *                        domain if necessary.
     * \return True if the particle can be successfully placed into the container,
     * false otherwise. */
    protected void put_locate_block(int[] ijk,double[] x,double[] y,double[] z) {
        // Remap particle in the z direction if necessary
        int k=step_int(z[0]*zsp);
        if(k<0||k>=nz) {
            int ak=step_div(k,nz);
            z[0]-=ak*bz;y[0]-=ak*byz;x[0]-=ak*bxz;k-=ak*nz;
        }

        // Remap particle in the y direction if necessary
        int j=step_int(y[0]*ysp);
        if(j<0||j>=ny) {
            int aj=step_div(j,ny);
            y[0]-=aj*by;x[0]-=aj*bxy;j-=aj*ny;
        }

        // Remap particle in the x direction if necessary
        ijk[0]=step_int(x[0]*xsp);
        if(ijk[0]<0||ijk[0]>=nx) {
            int ai=step_div(ijk[0],nx);
            x[0]-=ai*bx;ijk[0]-=ai*nx;
        }

        // Compute the block index and check memory allocation
        j+=ey;k+=ez;
        ijk[0]+=nx*(j+oy*k);
        if(co[ijk[0]]==mem[ijk[0]]) add_particle_memory(ijk[0]);
    }

    /** Takes a particle position vector and computes the region index into which
     * it should be stored. If the container is periodic, then the routine also
     * maps the particle position to ensure it is in the primary domain. If the
     * container is not periodic, the routine bails out.
     * \param[out] ijk the region index.
     * \param[in,out] (x,y,z) the particle position, remapped into the primary
     *                        domain if necessary.
     * \param[out] (ai,aj,ak) the periodic image displacement that the particle is
     *                        in, with (0,0,0) corresponding to the primary domain.
     * \return True if the particle can be successfully placed into the container,
     * false otherwise. */
    protected void put_locate_block(int[] ijk,double[] x,double[] y,double[] z,int[] ai,int[] aj,int[] ak) {
        // Remap particle in the z direction if necessary
        int k=step_int(z[0]*zsp);
        if(k<0||k>=nz) {
            ak[0]=step_div(k,nz);
            z[0]-=ak[0]*bz;y[0]-=ak[0]*byz;x[0]-=ak[0]*bxz;k-=ak[0]*nz;
        } else ak[0]=0;

        // Remap particle in the y direction if necessary
        int j=step_int(y[0]*ysp);
        if(j<0||j>=ny) {
            aj[0]=step_div(j,ny);
            y[0]-=aj[0]*by;x[0]-=aj[0]*bxy;j-=aj[0]*ny;
        } else aj[0]=0;

        // Remap particle in the x direction if necessary
        ijk[0]=step_int(x[0]*xsp);
        if(ijk[0]<0||ijk[0]>=nx) {
            ai[0]=step_div(ijk[0],nx);
            x[0]-=ai[0]*bx;ijk[0]-=ai[0]*nx;
        } else ai[0]=0;

        // Compute the block index and check memory allocation
        j+=ey;k+=ez;
        ijk[0]+=nx*(j+oy*k);
        if(co[ijk[0]]==mem[ijk[0]]) add_particle_memory(ijk[0]);
    }
    /** Creates particles within an image block by copying them
     * from the primary domain and shifting them. If the given
     * block is aligned with the primary domain in the z-direction,
     * the routine calls the simpler create_side_image routine
     * where the image block may comprise of particles from up to
     * two primary blocks. Otherwise is calls the more complex
     * create_vertical_image where the image block may comprise of
     * particles from up to four primary blocks.
     * \param[in] (di,dj,dk) the coordinates of the image block to
     *                       create. */
    protected void create_periodic_image(int di,int dj,int dk) {
        if(di<0||di>=nx||dj<0||dj>=oy||dk<0||dk>=oz)
            Common.voro_fatal_error("Constructing periodic image for nonexistent point",Config.Voropp.INTERNAL_ERROR);
        if(dk>=ez&&dk<wz) {
            if(dj<ey||dj>=wy) create_side_image(di,dj,dk);
        } else create_vertical_image(di,dj,dk);
    }

    /** Creates particles within an image block that is aligned with the primary
     * domain in the z axis. In this case, the image block may be comprised of
     * particles from two primary blocks. The routine considers these two primary
     * blocks, and adds the needed particles to the image. The remaining particles
     * from the primary blocks are also filled into the neighboring images.
     * \param[in] (di,dj,dk) the index of the block to consider. The z index must
     *			 satisfy ez<=dk<wz. */
    protected void create_side_image(int di,int dj,int dk) {
        int l,dijk=di+nx*(dj+oy*dk),odijk,ima=step_div(dj-ey,ny);
        int qua=di+step_int(-ima*bxy*xsp),quadiv=step_div(qua,nx);
        int fi=qua-quadiv*nx,fijk=fi+nx*(dj-ima*ny+oy*dk);
        double dis=ima*bxy+quadiv*bx,switchx=di*boxx-ima*bxy-quadiv*bx,adis;

        // Left image computation
        if((img[dijk]&1)==0) {
            if(di>0) {
                odijk=dijk-1;adis=dis;
            } else {
                odijk=dijk+nx-1;adis=dis+bx;
            }
            img[odijk]|=2;
            for(l=0;l<co[fijk];l++) {
                if(p[fijk][ps*l]>switchx) put_image(dijk,fijk,l,dis,by*ima,0);
                else put_image(odijk,fijk,l,adis,by*ima,0);
            }
        }

        // Right image computation
        if((img[dijk]&2)==0) {
            if(fi==nx-1) {
                fijk+=1-nx;switchx+=(1-nx)*boxx;dis+=bx;
            } else {
                fijk++;switchx+=boxx;
            }
            if(di==nx-1) {
                odijk=dijk-nx+1;adis=dis-bx;
            } else {
                odijk=dijk+1;adis=dis;
            }
            img[odijk]|=1;
            for(l=0;l<co[fijk];l++) {
                if(p[fijk][ps*l]<switchx) put_image(dijk,fijk,l,dis,by*ima,0);
                else put_image(odijk,fijk,l,adis,by*ima,0);
            }
        }

        // All contributions to the block now added, so set both two bits of
        // the image information
        img[dijk]=3;
    }

    /** Creates particles within an image block that is not aligned with the
     * primary domain in the z axis. In this case, the image block may be comprised
     * of particles from four primary blocks. The routine considers these four
     * primary blocks, and adds the needed particles to the image. The remaining
     * particles from the primary blocks are also filled into the neighboring
     * images.
     * \param[in] (di,dj,dk) the index of the block to consider. The z index must
     *			 satisfy dk<ez or dk>=wz. */
    protected void create_vertical_image(int di,int dj,int dk) {
        int l,dijk=di+nx*(dj+oy*dk),dijkl,dijkr,ima=step_div(dk-ez,nz);
        int qj=dj+step_int(-ima*byz*ysp),qjdiv=step_div(qj-ey,ny);
        int qi=di+step_int((-ima*bxz-qjdiv*bxy)*xsp),qidiv=step_div(qi,nx);
        int fi=qi-qidiv*nx,fj=qj-qjdiv*ny,fijk=fi+nx*(fj+oy*(dk-ima*nz)),fijk2;
        double disy=ima*byz+qjdiv*by,switchy=(dj-ey)*boxy-ima*byz-qjdiv*by;
        double disx=ima*bxz+qjdiv*bxy+qidiv*bx,switchx=di*boxx-ima*bxz-qjdiv*bxy-qidiv*bx;
        double switchx2,disxl,disxr,disx2,disxr2;

        if(di==0) {dijkl=dijk+nx-1;disxl=disx+bx;}
        else {dijkl=dijk-1;disxl=disx;}

        if(di==nx-1) {dijkr=dijk-nx+1;disxr=disx-bx;}
        else {dijkr=dijk+1;disxr=disx;}

        // Down-left image computation
        boolean y_exist=dj!=0;
        if((img[dijk]&1)==0) {
            img[dijkl]|=2;
            if(y_exist) {
                img[dijkl-nx]|=8;
                img[dijk-nx]|=4;
            }
            for(l=0;l<co[fijk];l++) {
                if(p[fijk][ps*l+1]>switchy) {
                    if(p[fijk][ps*l]>switchx) put_image(dijk,fijk,l,disx,disy,bz*ima);
                    else put_image(dijkl,fijk,l,disxl,disy,bz*ima);
                } else {
                    if(!y_exist) continue;
                    if(p[fijk][ps*l]>switchx) put_image(dijk-nx,fijk,l,disx,disy,bz*ima);
                    else put_image(dijkl-nx,fijk,l,disxl,disy,bz*ima);
                }
            }
        }

        // Down-right image computation
        if((img[dijk]&2)==0) {
            if(fi==nx-1) {
                fijk2=fijk+1-nx;switchx2=switchx+(1-nx)*boxx;disx2=disx+bx;disxr2=disxr+bx;
            } else {
                fijk2=fijk+1;switchx2=switchx+boxx;disx2=disx;disxr2=disxr;
            }
            img[dijkr]|=1;
            if(y_exist) {
                img[dijkr-nx]|=4;
                img[dijk-nx]|=8;
            }
            for(l=0;l<co[fijk2];l++) {
                if(p[fijk2][ps*l+1]>switchy) {
                    if(p[fijk2][ps*l]>switchx2) put_image(dijkr,fijk2,l,disxr2,disy,bz*ima);
                    else put_image(dijk,fijk2,l,disx2,disy,bz*ima);
                } else {
                    if(!y_exist) continue;
                    if(p[fijk2][ps*l]>switchx2) put_image(dijkr-nx,fijk2,l,disxr2,disy,bz*ima);
                    else put_image(dijk-nx,fijk2,l,disx2,disy,bz*ima);
                }
            }
        }

        // Recomputation of some intermediate quantities for boundary cases
        if(fj==wy-1) {
            fijk+=nx*(1-ny)-fi;
            switchy+=(1-ny)*boxy;
            disy+=by;
            qi=di+step_int(-(ima*bxz+(qjdiv+1)*bxy)*xsp);
            int dqidiv=step_div(qi,nx)-qidiv;qidiv+=dqidiv;
            fi=qi-qidiv*nx;
            fijk+=fi;
            disx+=bxy+bx*dqidiv;
            disxl+=bxy+bx*dqidiv;
            disxr+=bxy+bx*dqidiv;
            switchx-=bxy+bx*dqidiv;
        } else {
            fijk+=nx;switchy+=boxy;
        }

        // Up-left image computation
        y_exist=dj!=oy-1;
        if((img[dijk]&4)==0) {
            img[dijkl]|=8;
            if(y_exist) {
                img[dijkl+nx]|=2;
                img[dijk+nx]|=1;
            }
            for(l=0;l<co[fijk];l++) {
                if(p[fijk][ps*l+1]>switchy) {
                    if(!y_exist) continue;
                    if(p[fijk][ps*l]>switchx) put_image(dijk+nx,fijk,l,disx,disy,bz*ima);
                    else put_image(dijkl+nx,fijk,l,disxl,disy,bz*ima);
                } else {
                    if(p[fijk][ps*l]>switchx) put_image(dijk,fijk,l,disx,disy,bz*ima);
                    else put_image(dijkl,fijk,l,disxl,disy,bz*ima);
                }
            }
        }

        // Up-right image computation
        if((img[dijk]&8)==0) {
            if(fi==nx-1) {
                fijk2=fijk+1-nx;switchx2=switchx+(1-nx)*boxx;disx2=disx+bx;disxr2=disxr+bx;
            } else {
                fijk2=fijk+1;switchx2=switchx+boxx;disx2=disx;disxr2=disxr;
            }
            img[dijkr]|=4;
            if(y_exist) {
                img[dijkr+nx]|=1;
                img[dijk+nx]|=2;
            }
            for(l=0;l<co[fijk2];l++) {
                if(p[fijk2][ps*l+1]>switchy) {
                    if(!y_exist) continue;
                    if(p[fijk2][ps*l]>switchx2) put_image(dijkr+nx,fijk2,l,disxr2,disy,bz*ima);
                    else put_image(dijk+nx,fijk2,l,disx2,disy,bz*ima);
                } else {
                    if(p[fijk2][ps*l]>switchx2) put_image(dijkr,fijk2,l,disxr2,disy,bz*ima);
                    else put_image(dijk,fijk2,l,disx2,disy,bz*ima);
                }
            }
        }

        // All contributions to the block now added, so set all four bits of
        // the image information
        img[dijk]=15;
    }

    /** Copies a particle position from the primary domain into an image block.
     * \param[in] reg the block index within the primary domain that the particle
     *                is within.
     * \param[in] fijk the index of the image block.
     * \param[in] l the index of the particle entry within the primary block.
     * \param[in] (dx,dy,dz) the displacement vector to add to the particle. */
    protected void put_image(int reg,int fijk,int l,double dx,double dy,double dz) {
        if(co[reg]==mem[reg]) add_particle_memory(reg);
        int p1 = ps*co[reg]; // p[reg]
        int p2 = ps*l; // p[fijk]
        p[reg][p1] = p[fijk][p2] + dx;
        p[reg][p1+1] = p[fijk][p2+1] + dy;
        p[reg][p1+2] = p[fijk][p2+2] + dz;
        if(ps==4) p[reg][p1+3] = p[fijk][p2+3];
        id[reg][co[reg]++]=id[fijk][l];
    }

    /** Takes a position vector and remaps it into the primary domain.
     * \param[out] (ai,aj,ak) the periodic image displacement that the vector is in,
     *                        with (0,0,0) corresponding to the primary domain.
     * \param[out] (ci,cj,ck) the index of the block that the position vector is
     *                        within, once it has been remapped.
     * \param[in,out] (x,y,z) the position vector to consider, which is remapped
     *                        into the primary domain during the routine.
     * \param[out] ijk the block index that the vector is within. */
    protected  void remap(int[] ai,int[] aj,int[] ak,int[] ci,int[] cj,int[] ck,double[] x,double[] y,double[] z,int[] ijk) {
        // Remap particle in the z direction if necessary
        ck[0]=step_int(z[0]*zsp);
        if(ck[0]<0||ck[0]>=nz) {
            ak[0]=step_div(ck[0],nz);
            z[0]-=ak[0]*bz;
            y[0]-=ak[0]*byz;
            x[0]-=ak[0]*bxz;
            ck[0]-=ak[0]*nz;
        } else ak[0]=0;

        // Remap particle in the y direction if necessary
        cj[0]=step_int(y[0]*ysp);
        if(cj[0]<0||cj[0]>=ny) {
            aj[0]=step_div(cj[0],ny);
            y[0]-=aj[0]*by;
            x[0]-=aj[0]*bxy;
            cj[0]-=aj[0]*ny;
        } else aj[0]=0;

        // Remap particle in the x direction if necessary
        ci[0]=step_int(x[0]*xsp);
        if(ci[0]<0||ci[0]>=nx) {
            ai[0]=step_div(ci[0],nx);
            x[0]-=ai[0]*bx;
            ci[0]-=ai[0]*nx;
        } else ai[0]=0;

        cj[0]+=ey;ck[0]+=ez;
        ijk[0]=ci[0]+nx*(cj[0]+oy*ck[0]);
    }

}
