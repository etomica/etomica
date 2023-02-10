// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;

/** \brief Class for representing a particle system in a three-dimensional
 * rectangular box.
 *
 * This class represents a system of particles in a three-dimensional
 * rectangular box. Any combination of non-periodic and periodic coordinates
 * can be used in the three coordinate directions. The class is not intended
 * for direct use, but instead forms the base of the container and
 * container_poly classes that add specialized routines for computing the
 * regular and radical Voronoi tessellations respectively. It contains routines
 * that are commonly between these two classes, such as those for drawing the
 * domain, and placing particles within the internal data structure.
 *
 * The class is derived from the wall_list class, which encapsulates routines
 * for associating walls with the container, and the voro_base class, which
 * encapsulates routines about the underlying computational grid. */
public abstract class ContainerBase extends ContainerBaseBase implements Radius {

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
    /** This array holds the maximum amount of particle memory for
     * each computational box of the container. If the number of
     * particles in a particular box ever approaches this limit,
     * more is allocated using the add_particle_memory() function.
     */
    public int[] mem;

    public final ArrayList<Wall> wl = new ArrayList<>();

    /** The class constructor sets up the geometry of container, initializing the
     * minimum and maximum coordinates in each direction, and setting whether each
     * direction is periodic or not. It divides the container into a rectangular
     * grid of blocks, and allocates memory for each of these for storing particle
     * positions and IDs.
     * \param[in] (ax_,bx_) the minimum and maximum x coordinates.
     * \param[in] (ay_,by_) the minimum and maximum y coordinates.
     * \param[in] (az_,bz_) the minimum and maximum z coordinates.
     * \param[in] (nx_,ny_,nz_) the number of grid blocks in each of the three
     *			    coordinate directions.
     * \param[in] (xperiodic_,yperiodic_,zperiodic_) flags setting whether the
     *                                               container is periodic in each
     *                                               coordinate direction.
     * \param[in] init_mem the initial memory allocation for each block.
     * \param[in] ps_ the number of floating point entries to store for each
     *                particle. */
    public ContainerBase(double ax_,double bx_,double ay_,double by_,double az_,double bz_,
                int nx_,int ny_,int nz_,boolean xperiodic_,boolean yperiodic_,boolean zperiodic_,
                int init_mem,int ps_) {
        super(nx_,ny_,nz_,(bx_-ax_)/nx_,(by_-ay_)/ny_,(bz_-az_)/nz_, ps_);
        ax = ax_;
        bx = bx_;
        ay = ay_;
        by = by_;
        az = az_;
        bz = bz_;
        max_len_sq = (bx-ax)*(bx-ax)*(xperiodic_?0.25:1)+(by-ay)*(by-ay)*(yperiodic_?0.25:1)
                        +(bz-az)*(bz-az)*(zperiodic_?0.25:1);
        xperiodic = xperiodic_;
        yperiodic = yperiodic_;
        zperiodic = zperiodic_;
        id = new int[nxyz][];
        p  = new double[nxyz][];
        co = new int[nxyz];
        mem = new int[nxyz];

        int l;
        for(l=0;l<nxyz;l++) co[l]=0;
        for(l=0;l<nxyz;l++) mem[l]=init_mem;
        for(l=0;l<nxyz;l++) id[l]=new int[init_mem];
        for(l=0;l<nxyz;l++) p[l]=new double[ps*init_mem];
    }

    /** This function tests to see if a given vector lies within the container
     * bounds and any walls.
     * \param[in] (x,y,z) the position vector to be tested.
     * \return True if the point is inside the container, false if the point is
     *         outside. */
    public boolean point_inside(double x,double y,double z) {
        if(x<ax||x>bx||y<ay||y>by||z<az||z>bz) return false;
        return point_inside_walls(x,y,z);
    }

    /** Outputs the a list of all the container regions along with the number of
     * particles stored within each. */
    public void region_count() {
        for(int k=0,cop=0;k<nz;k++) for(int j=0;j<ny;j++) for(int i=0;i<nx;i++,cop++)
            System.out.printf("Region (%d,%d,%d): %d particles\n",i,j,k,co[cop]);
    }

    public boolean point_inside_walls(double x,double y,double z) {
        for (int i=0; i<wl.size(); i++) if (!wl.get(i).point_inside(x,y,z)) return false;
        return true;
    }

    /** Cuts a Voronoi cell by all of the walls currently on
     * the list.
     * \param[in] c a reference to the Voronoi cell class.
     * \param[in] (x,y,z) the position of the cell.
     * \return True if the cell still exists, false if the cell is
     * deleted. */
    public boolean apply_walls(VoronoiCellBase c, double x, double y, double z) {
        if (c instanceof VoronoiCell) {
            for (Wall w : wl) if (!(w.cut_cell((VoronoiCell) c, x, y, z))) return false;
        }
        else if (c instanceof VoronoiCellNeighbor) {
            for (Wall w : wl) if (!(w.cut_cell((VoronoiCellNeighbor) c, x, y, z))) return false;
        }
        return true;
    }

    /** Initializes the Voronoi cell prior to a compute_cell
     * operation for a specific particle being carried out by a
     * voro_compute class. The cell is initialized to fill the
     * entire container. For non-periodic coordinates, this is set
     * by the position of the walls. For periodic coordinates, the
     * space is equally divided in either direction from the
     * particle's initial position. Plane cuts made by any walls
     * that have been added are then applied to the cell.
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
    public boolean initialize_voronoicell(VoronoiCellBase[] c,int ijk,int q,int ci,int cj,int ck,
            int[] i,int[] j,int[] k,double[] x,double[] y,double[] z,int[] disp) {
        double x1,x2,y1,y2,z1,z2;
        x[0] = p[ijk][ps*q];
        y[0] = p[ijk][ps*q+1];
        z[0] = p[ijk][ps*q+2];
        if(xperiodic) {x1=-(x2=0.5*(bx-ax));i[0]=nx;} else {x1=ax-x[0];x2=bx-x[0];i[0]=ci;}
        if(yperiodic) {y1=-(y2=0.5*(by-ay));j[0]=ny;} else {y1=ay-y[0];y2=by-y[0];j[0]=cj;}
        if(zperiodic) {z1=-(z2=0.5*(bz-az));k[0]=nz;} else {z1=az-z[0];z2=bz-z[0];k[0]=ck;}
        c[0].init(x1,x2,y1,y2,z1,z2);
        if(!apply_walls(c[0],x[0],y[0],z[0])) return false;
        disp[0]=ijk-i[0]-nx*(j[0]+ny*k[0]);
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
     *		    find_voronoi_cell routine. */
    public void initialize_search(int ci,int cj,int ck,int ijk,int[] i,int[] j,int[] k,int[] disp) {
        i[0]=xperiodic?nx:ci;
        j[0]=yperiodic?ny:cj;
        k[0]=zperiodic?nz:ck;
        disp[0]=ijk-i[0]-nx*(j[0]+ny*k[0]);
    }
    /** Returns the position of a particle currently being computed
     * relative to the computational block that it is within. It is
     * used to select the optimal worklist entry to use.
     * \param[in] (x,y,z) the position of the particle.
     * \param[in] (ci,cj,ck) the block that the particle is within.
     * \param[out] (fx,fy,fz) the position relative to the block.
     */
    public void frac_pos(double x,double y,double z,double ci,double cj,double ck,
                    double[] fx,double[] fy,double[] fz) {
        fx[0]=x-ax-boxx*ci;
        fy[0]=y-ay-boxy*cj;
        fz[0]=z-az-boxz*ck;
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
     * 		    find_voronoi_cell and compute_cell routines.
     * \return The block index. */
    public int region_index(int ci,int cj,int ck,int ei,int ej,int ek,double[] qx,double[] qy,double[] qz,int[] disp) {
        if(xperiodic) {if(ci+ei<nx) {ei+=nx;qx[0]=-(bx-ax);} else if(ci+ei>=(nx<<1)) {ei-=nx;qx[0]=bx-ax;} else qx[0]=0;}
        if(yperiodic) {if(cj+ej<ny) {ej+=ny;qy[0]=-(by-ay);} else if(cj+ej>=(ny<<1)) {ej-=ny;qy[0]=by-ay;} else qy[0]=0;}
        if(zperiodic) {if(ck+ek<nz) {ek+=nz;qz[0]=-(bz-az);} else if(ck+ek>=(nz<<1)) {ek-=nz;qz[0]=bz-az;} else qz[0]=0;}
        return disp[0]+ei+nx*(ej+ny*ek);
    }

    /** Draws an outline of the domain in gnuplot format.
     * \param[in] fp the file handle to write to. */
    public void draw_domain_gnuplot() {
        draw_domain_gnuplot(System.out);
    }
    public void draw_domain_gnuplot(OutputStream fp) {
        try {
            fp.write(String.format("%g %g %g\n%g %g %g\n%g %g %g\n%g %g %g\n", ax, ay, az, bx, ay, az, bx, by, az, ax, by, az).getBytes());
            fp.write(String.format("%g %g %g\n%g %g %g\n%g %g %g\n%g %g %g\n", ax, by, bz, bx, by, bz, bx, ay, bz, ax, ay, bz).getBytes());
            fp.write(String.format("%g %g %g\n\n%g %g %g\n%g %g %g\n\n", ax, by, bz, ax, ay, az, ax, ay, bz).getBytes());
            fp.write(String.format("%g %g %g\n%g %g %g\n\n%g %g %g\n%g %g %g\n\n", bx, ay, az, bx, ay, bz, bx, by, az, bx, by, bz).getBytes());
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }
    /** Draws an outline of the domain in Gnuplot format.
     * \param[in] filename the filename to write to. */
    public void draw_domain_gnuplot(String filename) {
        try {
            FileOutputStream fp = new FileOutputStream(filename);
            draw_domain_gnuplot(fp);
            fp.close();
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }

    /** Draws an outline of the domain in POV-Ray format.
     * \param[in] fp the file handle to write to. */
    public void draw_domain_pov() {
        draw_domain_pov(System.out);
    }
    public void draw_domain_pov(OutputStream fp) {
        try {
            fp.write(String.format("cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n" +
                    "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n", ax, ay, az, bx, ay, az, ax, by, az, bx, by, az).getBytes());
            fp.write(String.format("cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n" +
                    "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n", ax, by, bz, bx, by, bz, ax, ay, bz, bx, ay, bz).getBytes());
            fp.write(String.format("cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n" +
                    "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n", ax, ay, az, ax, by, az, bx, ay, az, bx, by, az).getBytes());
            fp.write(String.format("cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n" +
                    "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n", bx, ay, bz, bx, by, bz, ax, ay, bz, ax, by, bz).getBytes());
            fp.write(String.format("cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n" +
                    "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n", ax, ay, az, ax, ay, bz, bx, ay, az, bx, ay, bz).getBytes());
            fp.write(String.format("cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n" +
                    "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n", bx, by, az, bx, by, bz, ax, by, az, ax, by, bz).getBytes());
            fp.write(String.format("sphere{<%g,%g,%g>,rr}\nsphere{<%g,%g,%g>,rr}\n" +
                    "sphere{<%g,%g,%g>,rr}\nsphere{<%g,%g,%g>,rr}\n", ax, ay, az, bx, ay, az, ax, by, az, bx, by, az).getBytes());
            fp.write(String.format("sphere{<%g,%g,%g>,rr}\nsphere{<%g,%g,%g>,rr}\n" +
                    "sphere{<%g,%g,%g>,rr}\nsphere{<%g,%g,%g>,rr}\n", ax, ay, bz, bx, ay, bz, ax, by, bz, bx, by, bz).getBytes());
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }
    /** Draws an outline of the domain in Gnuplot format.
     * \param[in] filename the filename to write to. */
    public void draw_domain_pov(String filename) {
        try {
            OutputStream fp = new FileOutputStream(filename);
            draw_domain_pov(fp);
            fp.close();
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }
    /** Sums up the total number of stored particles.
     * \return The number of particles. */
    public int total_particles() {
        int tp=co[0];
        for(int cop=1;cop<nxyz;cop++) tp+=co[cop];
        return tp;
    }

    /** Increase memory for a particular region.
     * \param[in] i the index of the region to reallocate. */
    protected void add_particle_memory(int i) {
        int l,nmem=mem[i]<<1;

        // Carry out a check on the memory allocation size, and
        // print a status message if requested
        if(nmem>Config.max_particle_memory)
            Common.voro_fatal_error("Absolute maximum memory allocation exceeded",Config.Voropp.MEMORY_ERROR);
        if (Config.VOROPP_VERBOSE >= 3) {
            System.err.printf("Particle memory in region %d scaled up to %d\n", i, nmem);
        }

        // Allocate new memory and copy in the contents of the old arrays
        int[] idp=new int[nmem];
        for(l=0;l<co[i];l++) idp[l]=id[i][l];
        double[] pp=new double[ps*nmem];
        for(l=0;l<ps*co[i];l++) pp[l]=p[i][l];

        // Update pointers and delete old arrays
        mem[i]=nmem;
        id[i]=idp;
        p[i]=pp;
    }

    /** This routine takes a particle position vector, tries to remap it into the
     * primary domain. If successful, it computes the region into which it can be
     * stored and checks that there is enough memory within this region to store
     * it.
     * \param[out] ijk the region index.
     * \param[in,out] (x,y,z) the particle position, remapped into the primary
     *                        domain if necessary.
     * \return True if the particle can be successfully placed into the container,
     * false otherwise. */
    protected boolean put_locate_block(int[] ijk,double[] x,double[] y,double[] z) {
        if(put_remap(ijk,x,y,z)) {
            if(co[ijk[0]]==mem[ijk[0]]) add_particle_memory(ijk[0]);
            return true;
        }
        if (Config.VOROPP_REPORT_OUT_OF_BOUNDS) {
            System.err.printf("Out of bounds: (x,y,z)=(%g,%g,%g)\n", x[0], y[0], z[0]);
        }
        return false;
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
    protected boolean put_remap(int[] ijk,double[] x,double[] y,double[] z) {
        int l;
        ijk[0]=step_int((x[0]-ax)*xsp);
        if(xperiodic) {l=step_mod(ijk[0],nx);x[0]+=boxx*(l-ijk[0]);ijk[0]=l;}
        else if(ijk[0]<0||ijk[0]>=nx) return false;

        int j=step_int((y[0]-ay)*ysp);
        if(yperiodic) {l=step_mod(j,ny);y[0]+=boxy*(l-j);j=l;}
        else if(j<0||j>=ny) return false;

        int k=step_int((z[0]-az)*zsp);
        if(zperiodic) {l=step_mod(k,nz);z[0]+=boxz*(l-k);k=l;}
        else if(k<0||k>=nz) return false;

        ijk[0]+=nx*j+nxy*k;
        return true;
    }

    /** Takes a position vector and attempts to remap it into the primary domain.
     * \param[out] (ai,aj,ak) the periodic image displacement that the vector is in,
     *                       with (0,0,0) corresponding to the primary domain.
     * \param[out] (ci,cj,ck) the index of the block that the position vector is
     *                        within, once it has been remapped.
     * \param[in,out] (x,y,z) the position vector to consider, which is remapped
     *                        into the primary domain during the routine.
     * \param[out] ijk the block index that the vector is within.
     * \return True if the particle is within the container or can be remapped into
     * it, false if it lies outside of the container bounds. */
    protected boolean remap(int[] ai,int[] aj,int[] ak,int[] ci,int[] cj,int[] ck,double[] x,double[] y,double[] z,int[] ijk) {
        ci[0]=step_int((x[0]-ax)*xsp);
        if(ci[0]<0||ci[0]>=nx) {
            if(xperiodic) {ai[0]=step_div(ci[0],nx);x[0]-=ai[0]*(bx-ax);ci[0]-=ai[0]*nx;}
            else return false;
        } else ai[0]=0;

        cj[0]=step_int((y[0]-ay)*ysp);
        if(cj[0]<0||cj[0]>=ny) {
            if(yperiodic) {aj[0]=step_div(cj[0],ny);y[0]-=aj[0]*(by-ay);cj[0]-=aj[0]*ny;}
            else return false;
        } else aj[0]=0;

        ck[0]=step_int((z[0]-az)*zsp);
        if(ck[0]<0||ck[0]>=nz) {
            if(zperiodic) {ak[0]=step_div(ck[0],nz);z[0]-=ak[0]*(bz-az);ck[0]-=ak[0]*nz;}
            else return false;
        } else ak[0]=0;

        ijk[0]=ci[0]+nx*cj[0]+nxy*ck[0];
        return true;
    }
}
