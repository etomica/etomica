/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.util.voro;

import java.io.*;
import java.util.Arrays;

/** \brief Extension of the container_base class for computing regular Voronoi
 * tessellations.
 *
 * This class is an extension of the container_base class that has routines
 * specifically for computing the regular Voronoi tessellation with no
 * dependence on particle radii. */
public class Container extends ContainerBase {

    /** The class constructor sets up the geometry of container.
     * \param[in] (ax_,bx_) the minimum and maximum x coordinates.
     * \param[in] (ay_,by_) the minimum and maximum y coordinates.
     * \param[in] (az_,bz_) the minimum and maximum z coordinates.
     * \param[in] (nx_,ny_,nz_) the number of grid blocks in each of the three
     *                       coordinate directions.
     * \param[in] (xperiodic_,yperiodic_,zperiodic_) flags setting whether the
     *                                               container is periodic in each
     *                                               coordinate direction.
     * \param[in] init_mem the initial memory allocation for each block. */
    public Container(double ax_,double bx_,double ay_,double by_,double az_,double bz_,
              int nx_,int ny_,int nz_,boolean xperiodic_,boolean yperiodic_,boolean zperiodic_,int init_mem) {
	    super(ax_,bx_,ay_,by_,az_,bz_,nx_,ny_,nz_,xperiodic_,yperiodic_,zperiodic_,init_mem,3);
        radius = new RadiusMono();
        vc = new VoroCompute(this,xperiodic_?2*nx_+1:nx_,yperiodic_?2*ny_+1:ny_,zperiodic_?2*nz_+1:nz_);
    }

    /** Clears a container of particles. */
    public void clear() {
        Arrays.fill(co, 0);
    }

    /** Put a particle into the correct region of the container.
     * \param[in] n the numerical ID of the inserted particle.
     * \param[in] (x,y,z) the position vector of the inserted particle. */
    public void put(int n,double x,double y,double z) {
        int[] ijk = new int[1];
        double[] xout = new double[]{x};
        double[] yout = new double[]{y};
        double[] zout = new double[]{z};
        if(put_locate_block(ijk,xout,yout,zout)) {
            id[ijk[0]][co[ijk[0]]]=n;
            double[] pijk = p[ijk[0]];
            int pp = 3*co[ijk[0]]; // pijk
            co[ijk[0]]++;
            pijk[pp+0] = xout[0];
            pijk[pp+1] = yout[0];
            pijk[pp+2] = zout[0];
        }
    }

    /** Put a particle into the correct region of the container, also recording
     * into which region it was stored.
     * \param[in] vo the ordering class in which to record the region.
     * \param[in] n the numerical ID of the inserted particle.
     * \param[in] (x,y,z) the position vector of the inserted particle. */
    public void put(ParticleOrder vo,int n,double x,double y,double z) {
        int[] ijk = new int[1];
        double[] xout = new double[]{x};
        double[] yout = new double[]{y};
        double[] zout = new double[]{z};
        if(put_locate_block(ijk,xout,yout,zout)) {
            id[ijk[0]][co[ijk[0]]]=n;
            vo.add(ijk[0],co[ijk[0]]);
            double[] pijk = p[ijk[0]];
            int pp = 3*co[ijk[0]];
            co[ijk[0]]++;
            pijk[pp+0] = xout[0];
            pijk[pp+1] = yout[0];
            pijk[pp+2] = zout[0];
        }
    }

    /** Import a list of particles from an open file stream into the container.
     * Entries of four numbers (Particle ID, x position, y position, z position)
     * are searched for. If the file cannot be successfully read, then the routine
     * causes a fatal error.
     * \param[in] fp the file handle to read from. */
    public void import_() {
        import_(System.in);
    }
    public void import_(InputStream fp) {
        BufferedReader bufReader = new BufferedReader(new InputStreamReader(fp));

        try {
            String line = null;
            while ((line = bufReader.readLine()) != null) {
                String[] bits = line.trim().split("[ \t]+");
                if (bits.length != 4) Common.voro_fatal_error("File import error", Config.Voropp.FILE_ERROR);

                int i = Integer.parseInt(bits[0]);
                double x = Double.parseDouble(bits[1]);
                double y = Double.parseDouble(bits[2]);
                double z = Double.parseDouble(bits[3]);
                put(i, x, y, z);
            }
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }
    public void import_(ParticleOrder vo) {
        import_(vo, System.in);
    }
    public void import_(ParticleOrder vo,InputStream fp) {
        BufferedReader bufReader = new BufferedReader(new InputStreamReader(fp));

        try {
            String line = null;
            while ((line = bufReader.readLine()) != null) {
                String[] bits = line.trim().split("[ \t]+");
                if (bits.length != 4) Common.voro_fatal_error("File impor error", Config.Voropp.FILE_ERROR);

                int i = Integer.parseInt(bits[0]);
                double x = Double.parseDouble(bits[1]);
                double y = Double.parseDouble(bits[2]);
                double z = Double.parseDouble(bits[3]);
                put(vo, i, x, y, z);
            }
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }
    /** Imports a list of particles from an open file stream into
     * the container. Entries of four numbers (Particle ID, x
     * position, y position, z position) are searched for. If the
     * file cannot be successfully read, then the routine causes a
     * fatal error.
     * \param[in] filename the name of the file to open and read
     *                     from. */
    public void import_(String filename) {
        try {
            FileInputStream fp = new FileInputStream(filename);
            import_(fp);
            fp.close();
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }
    /** Imports a list of particles from an open file stream into
     * the container. Entries of four numbers (Particle ID, x
     * position, y position, z position) are searched for. In
     * addition, the order in which particles are read is saved
     * into an ordering class. If the file cannot be successfully
     * read, then the routine causes a fatal error.
     * \param[in,out] vo the ordering class to use.
     * \param[in] filename the name of the file to open and read
     *                     from. */
    public void import_(ParticleOrder vo,String filename) {
        try {
            FileInputStream fp = new FileInputStream(filename);
            import_(vo,fp);
            fp.close();
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }

    /** Computes all of the Voronoi cells in the container, but does nothing
     * with the output. It is useful for measuring the pure computation time
     * of the Voronoi algorithm, without any additional calculations such as
     * volume evaluation or cell output. */
    public void compute_all_cells() {
        CLoopAll vl = new CLoopAll(this);
        if(vl.start()) {
            VoronoiCell c = new VoronoiCell(this);
            VoronoiCell[] cout = new VoronoiCell[]{c};
            do {
                compute_cell(cout,vl);
            }
            while(vl.inc());
        }
    }

    /** Calculates all of the Voronoi cells and sums their volumes. In most cases
     * without walls, the sum of the Voronoi cell volumes should equal the volume
     * of the container to numerical precision.
     * \return The sum of all of the computed Voronoi volumes. */
    public double sum_cell_volumes() {
        VoronoiCell c = new VoronoiCell(this);
        VoronoiCell[] cout = new VoronoiCell[]{c};
        double vol=0;
        CLoopAll vl = new CLoopAll(this);
        if(vl.start()) {
            do {
                if(compute_cell(cout,vl)) vol+=cout[0].volume();
            }while(vl.inc());
        }
        return vol;
    }
    /** Dumps particle IDs and positions to a file.
     * \param[in] vl the loop class to use.
     * \param[in] fp a file handle to write to. */
    public void draw_particles(CLoopBase vl,OutputStream fp) {
        try {
            if (vl.start()) do {
                int pp = 3 * vl.q; //p[vl.ijk]
                fp.write(String.format("%d %g %g %g\n", id[vl.ijk][vl.q], p[vl.ijk][pp], p[vl.ijk][pp + 1], p[vl.ijk][pp + 2]).getBytes());
            } while (vl.inc());
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }
    /** Dumps all of the particle IDs and positions to a file.
     * \param[in] fp a file handle to write to. */
    public void draw_particles() {
        draw_particles(System.out);
    }
    public void draw_particles(OutputStream fp) {
        CLoopAll vl = new CLoopAll(this);
        draw_particles(vl,fp);
    }
    /** Dumps all of the particle IDs and positions to a file.
     * \param[in] filename the name of the file to write to. */
    public void draw_particles(String filename) {
        try {
            FileOutputStream fp = new FileOutputStream(filename);
            draw_particles(fp);
            fp.close();
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }
    /** Dumps particle positions in POV-Ray format.
     * \param[in] vl the loop class to use.
     * \param[in] fp a file handle to write to. */
    void draw_particles_pov(CLoopBase vl, OutputStream fp) {
        try {
            if (vl.start()) do {
                int pp = 3 * vl.q; // p[vl.ijk]
                fp.write(String.format("// id %d\nsphere{<%g,%g,%g>,s}\n",
                        id[vl.ijk][vl.q], p[vl.ijk][pp], p[vl.ijk][pp + 1], p[vl.ijk][pp + 2]).getBytes());
            } while (vl.inc());
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }
    /** Dumps all particle positions in POV-Ray format.
     * \param[in] fp a file handle to write to. */
    public void draw_particles_pov() {
        draw_particles_pov(System.out);
    }
    public void draw_particles_pov(OutputStream fp) {
        CLoopAll vl = new CLoopAll(this);
        draw_particles_pov(vl,fp);
    }
    /** Dumps all particle positions in POV-Ray format.
     * \param[in] filename the name of the file to write to. */
    public void draw_particles_pov(String filename) {
        try {
            FileOutputStream fp = new FileOutputStream(filename);
            draw_particles_pov(fp);
            fp.close();
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }
    /** Computes Voronoi cells and saves the output in gnuplot
     * format.
     * \param[in] vl the loop class to use.
     * \param[in] fp a file handle to write to. */
    void draw_cells_gnuplot(CLoopBase vl,OutputStream fp) {
        if(vl.start()) {
            VoronoiCell c = new VoronoiCell(this);
            VoronoiCell[] cout = new VoronoiCell[]{c};
            do {
                if(compute_cell(cout,vl)) {
                    int pp = ps*vl.q; // p[vl.ijk]
                    cout[0].draw_gnuplot(p[vl.ijk][pp],p[vl.ijk][pp+1],p[vl.ijk][pp+2],fp);
                }
            } while(vl.inc());
        }
    }
    /** Computes all Voronoi cells and saves the output in gnuplot
     * format.
     * \param[in] fp a file handle to write to. */
    public void draw_cells_gnuplot() {
        draw_cells_gnuplot(System.out);
    }
    public void draw_cells_gnuplot(OutputStream fp) {
        CLoopAll vl = new CLoopAll(this);
        draw_cells_gnuplot(vl,fp);
    }
    /** Computes all Voronoi cells and saves the output in gnuplot
     * format.
     * \param[in] filename the name of the file to write to. */
    public void draw_cells_gnuplot(String filename) {
        try {
            FileOutputStream fp = new FileOutputStream(filename);
            draw_cells_gnuplot(fp);
            fp.close();
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }
    /** Computes Voronoi cells and saves the output in POV-Ray
     * format.
     * \param[in] vl the loop class to use.
     * \param[in] fp a file handle to write to. */
    void draw_cells_pov(CLoopBase vl,OutputStream fp) {
        try {
            if (vl.start()) {
                VoronoiCell c = new VoronoiCell(this);
                VoronoiCell[] cout = new VoronoiCell[]{c};
                do {
                    if (compute_cell(cout, vl)) {
                        fp.write(String.format("// cell %d\n", id[vl.ijk][vl.q]).getBytes());
                        int pp = ps * vl.q; // p[vl.ijk]
                        cout[0].draw_pov(p[vl.ijk][pp], p[vl.ijk][pp + 1], p[vl.ijk][pp + 2], fp);
                    }
                } while (vl.inc());
            }
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }
    /** Computes all Voronoi cells and saves the output in POV-Ray
     * format.
     * \param[in] fp a file handle to write to. */
    public void draw_cells_pov() {
        draw_cells_pov(System.out);
    }
    public void draw_cells_pov(OutputStream fp) {
        CLoopAll vl = new CLoopAll(this);
        draw_cells_pov(vl,fp);
    }
    /** Computes all Voronoi cells and saves the output in POV-Ray
     * format.
     * \param[in] filename the name of the file to write to. */
    public void draw_cells_pov(String filename) {
        try {
            FileOutputStream fp = new FileOutputStream(filename);
            draw_cells_pov(fp);
            fp.close();
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }
    /** Computes the Voronoi cells and saves customized information
     * about them.
     * \param[in] vl the loop class to use.
     * \param[in] format the custom output string to use.
     * \param[in] fp a file handle to write to. */
    public void print_custom(CLoopBase vl,String format,OutputStream fp) {
        int ijk,q;
        if(contains_neighbor(format)) {
            if(vl.start()) {
                VoronoiCellNeighbor c = new VoronoiCellNeighbor(this);
                VoronoiCellNeighbor[] cout = new VoronoiCellNeighbor[]{c};
                do {
                    if(compute_cell(cout,vl)) {
                        ijk=vl.ijk;
                        q=vl.q;
                        int pp = ps*q; // p[ijk]
                        cout[0].output_custom(format,id[ijk][q],p[ijk][pp],p[ijk][pp+1],p[ijk][pp+2],Config.default_radius,fp);
                    }
                } while(vl.inc());
            }
        } else {
            if(vl.start()) {
                VoronoiCell c = new VoronoiCell(this);
                VoronoiCell[] cout = new VoronoiCell[]{c};
                do {
                    if(compute_cell(cout,vl)) {
                        ijk=vl.ijk;q=vl.q;
                        int pp = ps*q; //p[ijk]
                        cout[0].output_custom(format,id[ijk][q],p[ijk][pp],p[ijk][pp+1],p[ijk][pp+2],Config.default_radius,fp);
                    }
                } while(vl.inc());
            }
        }
    }

    /** Computes all the Voronoi cells and saves customized information about them.
     * \param[in] format the custom output string to use.
     * \param[in] fp a file handle to write to. */
    public void print_custom(String format) {
        print_custom(format, System.out);
    }
    void print_custom(String format, OutputStream fp) {
        CLoopAll vl = new CLoopAll(this);
        print_custom(vl,format,fp);
    }

    /** Computes all the Voronoi cells and saves customized information about them.
     * \param[in] format the custom output string to use.
     * \param[in] filename the name of the file to write to. */
    void print_custom(String format, String filename) {
        try {
            FileOutputStream fp = new FileOutputStream(filename);
            print_custom(format, fp);
            fp.close();
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }

    /** Takes a vector and finds the particle whose Voronoi cell contains that
     * vector. This is equivalent to finding the particle which is nearest to the
     * vector. Additional wall classes are not considered by this routine.
     * \param[in] (x,y,z) the vector to test.
     * \param[out] (rx,ry,rz) the position of the particle whose Voronoi cell
     *                        contains the vector. If the container is periodic,
     *                        this may point to a particle in a periodic image of
     *                        the primary domain.
     * \param[out] pid the ID of the particle.
     * \return True if a particle was found. If the container has no particles,
     * then the search will not find a Voronoi cell and false is returned. */
    public boolean find_voronoi_cell(double x,double y,double z,double[] rx,double[] ry,double[] rz,int[] pid) {
        int[] ai = new int[1], aj = new int[1], ak = new int[1], ci = new int[1], cj = new int[1], ck = new int[1], ijk = new int[1];
        ParticleRecord w = new ParticleRecord();
        double[] mrs = new double[1];

        // If the given vector lies outside the domain, but the container
        // is periodic, then remap it back into the domain
        double[] xout = new double[]{x}, yout = new double[]{y}, zout = new double[]{z};
        if(!remap(ai,aj,ak,ci,cj,ck,xout,yout,zout,ijk)) return false;
        x = xout[0];
        y = yout[0];
        z = zout[0];
        vc.find_voronoi_cell(x,y,z,ci[0],cj[0],ck[0],ijk[0],w,mrs);

        if(w.ijk!=-1) {

            // Assemble the position vector of the particle to be returned,
            // applying a periodic remapping if necessary
            if(xperiodic) {ci[0]+=w.di;if(ci[0]<0||ci[0]>=nx) ai[0]+=step_div(ci[0],nx);}
            if(yperiodic) {cj[0]+=w.dj;if(cj[0]<0||cj[0]>=ny) aj[0]+=step_div(cj[0],ny);}
            if(zperiodic) {ck[0]+=w.dk;if(ck[0]<0||ck[0]>=nz) ak[0]+=step_div(ck[0],nz);}
            rx[0]=p[w.ijk][3*w.l]+ai[0]*(bx-ax);
            ry[0]=p[w.ijk][3*w.l+1]+aj[0]*(by-ay);
            rz[0]=p[w.ijk][3*w.l+2]+ak[0]*(bz-az);
            pid[0]=id[w.ijk][w.l];
            return true;
        }

        // If no particle if found then just return false
        return false;
    }
    /** Computes the Voronoi cell for a particle currently being
     * referenced by a loop class.
     * \param[out] c a Voronoi cell class in which to store the
     * 		 computed cell.
     * \param[in] vl the loop class to use.
     * \return True if the cell was computed. If the cell cannot be
     * computed, if it is removed entirely by a wall or boundary
     * condition, then the routine returns false. */
    public boolean compute_cell(VoronoiCellBase[] c,CLoopBase vl) {
        return vc.compute_cell(c,vl.ijk,vl.q,vl.i,vl.j,vl.k);
    }
    /** Computes the Voronoi cell for given particle.
     * \param[out] c a Voronoi cell class in which to store the
     * 		 computed cell.
     * \param[in] ijk the block that the particle is within.
     * \param[in] q the index of the particle within the block.
     * \return True if the cell was computed. If the cell cannot be
     * computed, if it is removed entirely by a wall or boundary
     * condition, then the routine returns false. */
    public boolean compute_cell(VoronoiCellBase[] c,int ijk,int q) {
        int k=ijk/nxy,ijkt=ijk-nxy*k,j=ijkt/nx,i=ijkt-j*nx;
        return vc.compute_cell(c,ijk,q,i,j,k);
    }
    /** Computes the Voronoi cell for a ghost particle at a given
     * location.
     * \param[out] c a Voronoi cell class in which to store the
     * 		 computed cell.
     * \param[in] (x,y,z) the location of the ghost particle.
     * \return True if the cell was computed. If the cell cannot be
     * computed, if it is removed entirely by a wall or boundary
     * condition, then the routine returns false. */
    public boolean compute_ghost_cell(VoronoiCellBase[] c,double x,double y,double z) {
        int[] ijk = new int[1];
        double[] xout = new double[]{x}, yout = new double[]{y}, zout = new double[]{z};
        boolean plb = put_locate_block(ijk,xout,yout,zout);
        x = xout[0];
        y = yout[0];
        z = zout[0];
        if (plb) {
            int pp = 3*co[ijk[0]]; // p[ijk]
            co[ijk[0]]++;
            p[ijk[0]][pp] = x;
            p[ijk[0]][pp+1] = y;
            p[ijk[0]][pp+2] = z;
            boolean q=compute_cell(c,ijk[0],co[ijk[0]]-1);
            co[ijk[0]]--;
            return q;
        }
        return false;
    }
    private final VoroCompute vc;

}
