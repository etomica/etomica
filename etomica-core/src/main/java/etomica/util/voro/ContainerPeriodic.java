// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro;

import java.io.*;

/** \brief Extension of the container_periodic_base class for computing regular
 * Voronoi tessellations.
 *
 * This class is an extension of the container_periodic_base that has routines
 * specifically for computing the regular Voronoi tessellation with no
 * dependence on particle radii. */
public class ContainerPeriodic extends ContainerPeriodicBase {

    /** The class constructor sets up the geometry of container.
     * \param[in] (bx_) The x coordinate of the first unit vector.
     * \param[in] (bxy_,by_) The x and y coordinates of the second unit vector.
     * \param[in] (bxz_,byz_,bz_) The x, y, and z coordinates of the third unit
     *                            vector.
     * \param[in] (nx_,ny_,nz_) the number of grid blocks in each of the three
     *			    coordinate directions.
     * \param[in] init_mem_ the initial memory allocation for each block. */
    public ContainerPeriodic(double bx_,double bxy_,double by_,double bxz_,double byz_,double bz_,
                       int nx_,int ny_,int nz_,int init_mem_) {
       	super(bx_,bxy_,by_,bxz_,byz_,bz_,nx_,ny_,nz_,init_mem_,3);
        vc = new VoroCompute(this,2*nx_+1,2*ey+1,2*ez+1);
        radius = new RadiusMono();
    }

    /** Clears a container of particles. */
    public void clear() {
        for (int cop=0;cop<oxyz;cop++) co[cop]=0;
        for (int cp=0; cp<oxyz; cp++) img[cp] = 0;
    }

    /** Put a particle into the correct region of the container.
     * \param[in] n the numerical ID of the inserted particle.
     * \param[in] (x,y,z) the position vector of the inserted particle. */
    public void put(int n,double x,double y,double z) {
        int[] ijkout = new int[1];
        double[] xout = new double[]{x}, yout = new double[]{y}, zout = new double[]{z};
        put_locate_block(ijkout,xout,yout,zout);
        int ijk = ijkout[0];
        x = xout[0];
        y = yout[0];
        z = zout[0];
        for(int l=0;l<co[ijk];l++) Common.check_duplicate(n,x,y,z,id[ijk][l],p[ijk], 3*l);
        id[ijk][co[ijk]]=n;
        p[ijk][3*co[ijk]] = x;
        p[ijk][3*co[ijk]+1] = y;
        p[ijk][3*co[ijk]+2] = z;
        co[ijk]++;
    }

    /** Put a particle into the correct region of the container.
     * \param[in] n the numerical ID of the inserted particle.
     * \param[in] (x,y,z) the position vector of the inserted particle.
     * \param[out] (ai,aj,ak) the periodic image displacement that the particle is
     * 			  in, with (0,0,0) corresponding to the primary domain.
     */
    public void put(int n,double x,double y,double z,int[] ai,int[] aj,int[] ak) {
        int[] ijkout = new int[1];
        double[] xout = new double[]{x};
        double[] yout = new double[]{y};
        double[] zout = new double[]{z};
        put_locate_block(ijkout,xout,yout,zout,ai,aj,ak);
        int ijk = ijkout[0];
        x = xout[0];
        y = yout[0];
        z = zout[0];
        for(int l=0;l<co[ijk];l++) Common.check_duplicate(n,x,y,z,id[ijk][l],p[ijk], 3*l);
        id[ijk][co[ijk]]=n;
        p[ijk][3*co[ijk]] = x;
        p[ijk][3*co[ijk]+1] = y;
        p[ijk][3*co[ijk]+2] = z;
        co[ijk]++;
    }

    /** Put a particle into the correct region of the container, also recording
     * into which region it was stored.
     * \param[in] vo the ordering class in which to record the region.
     * \param[in] n the numerical ID of the inserted particle.
     * \param[in] (x,y,z) the position vector of the inserted particle. */
    public void put(ParticleOrder vo,int n,double x,double y,double z) {
        int[] ijkout = new int[1];
        double[] xout = new double[]{x};
        double[] yout = new double[]{y};
        double[] zout = new double[]{z};
        put_locate_block(ijkout,xout,yout,zout);
        int ijk = ijkout[0];
        x = xout[0];
        y = yout[0];
        z = zout[0];
        id[ijk][co[ijk]]=n;
        vo.add(ijk,co[ijk]);
        p[ijk][3*co[ijk]] = x;
        p[ijk][3*co[ijk]+1] = y;
        p[ijk][3*co[ijk]+2] = z;
        co[ijk]++;
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
        import_(null, fp);
    }

    /** Import a list of particles from an open file stream, also storing the order
     * of that the particles are read. Entries of four numbers (Particle ID, x
     * position, y position, z position) are searched for. If the file cannot be
     * successfully read, then the routine causes a fatal error.
     * \param[in,out] vo a reference to an ordering class to use.
     * \param[in] fp the file handle to read from. */
    public void import_(ParticleOrder vo) {
        import_(vo, System.in);
    }
    public void import_(ParticleOrder vo, InputStream fp) {
        BufferedReader bufReader = new BufferedReader(new InputStreamReader(fp));
        String line = null;
        try {
            while ((line = bufReader.readLine()) != null) {
                String[] bits = line.trim().split("[ \t]+");
                if (bits.length != 4) Common.voro_fatal_error("File impor error", Config.Voropp.FILE_ERROR);
                int i = Integer.parseInt(bits[0]);
                double x = Double.parseDouble(bits[1]);
                double y = Double.parseDouble(bits[2]);
                double z = Double.parseDouble(bits[3]);
                if (vo != null) put(vo, i, x, y, z);
                else put(i, x, y, z);
            }
        }
        catch(IOException ex) {
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
            import_(vo, fp);
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
        CLoopAllPeriodic vl = new CLoopAllPeriodic(this);
        if(vl.start()) {
            VoronoiCell c = new VoronoiCell(this);
            do {
                compute_cell(c,vl);
            }
            while(vl.inc());
        }
    }

    /** Calculates all of the Voronoi cells and sums their volumes. In most cases
     * without walls, the sum of the Voronoi cell volumes should equal the volume
     * of the container to numerical precision.
     * \return The sum of all of the computed Voronoi volumes. */
    public double sum_cell_volumes() {
        double vol=0;
        CLoopAllPeriodic vl = new CLoopAllPeriodic(this);
        if(vl.start()) {
            VoronoiCell c = new VoronoiCell(this);
            do {
                if(compute_cell(c,vl)) vol+=c.volume();
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
                double[] pijk = p[vl.ijk];
                fp.write(String.format("%d %g %g %g\n", id[vl.ijk][vl.q], pijk[pp], pijk[pp + 1], pijk[pp + 2]).getBytes());
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
        CLoopAllPeriodic vl = new CLoopAllPeriodic(this);
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
    public void draw_particles_pov(CLoopBase vl,OutputStream fp) {
        try {
            if (vl.start()) do {
                double[] pijk = p[vl.ijk];
                int pp = 3 * vl.q;
                fp.write(String.format("// id %d\nsphere{<%g,%g,%g>,s}\n",
                        id[vl.ijk][vl.q], pijk[pp], pijk[pp + 1], pijk[pp + 2]).getBytes());
            } while (vl.inc());
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }
    /** Dumps all particle positions in POV-Ray format.
     * \param[in] fp a file handle to write to. */
    public void draw_partcles_pov() {
        draw_particles_pov(System.out);
    }
    public void draw_particles_pov(OutputStream fp) {
        CLoopAllPeriodic vl = new CLoopAllPeriodic(this);
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
    void draw_cells_gnuplot(CLoopBase vl, OutputStream fp) {
        if(vl.start()) {
            do {
                VoronoiCell c = new VoronoiCell(this);
                VoronoiCell[] cout = new VoronoiCell[]{c};
                if(compute_cell(c,vl)) {
                    double[] pijk = p[vl.ijk];
                    int pp=+ps*vl.q;
                    cout[0].draw_gnuplot(pijk[pp],pijk[pp+1],pijk[pp+2],fp);
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
        CLoopAllPeriodic vl = new CLoopAllPeriodic(this);
        draw_cells_gnuplot(vl,fp);
    }
    /** Compute all Voronoi cells and saves the output in gnuplot
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
                    if (compute_cell(c, vl)) {
                        fp.write(String.format("// cell %d\n", id[vl.ijk][vl.q]).getBytes());
                        double[] pijk = p[vl.ijk];
                        int pp = ps * vl.q;
                        cout[0].draw_pov(pijk[pp], pijk[pp + 1], pijk[pp + 2], fp);
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
        CLoopAllPeriodic vl = new CLoopAllPeriodic(this);
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
    void print_custom(CLoopBase vl,String format,OutputStream fp) {
        int ijk,q;
        if(contains_neighbor(format)) {
            if(vl.start()) {
                VoronoiCellNeighbor c = new VoronoiCellNeighbor(this);
                do {
                    if(compute_cell(c,vl)) {
                        ijk=vl.ijk;
                        q=vl.q;
                        int pp=ps*q; //p[ijk]
                        c.output_custom(format,id[ijk][q],p[ijk][pp],p[ijk][pp+1],p[ijk][pp+2],Config.default_radius,fp);
                    }
                } while(vl.inc());
            }
        } else {
            if(vl.start()) {
                VoronoiCell c = new VoronoiCell(this);
                do {
                    if(compute_cell(c,vl)) {
                        ijk=vl.ijk;
                        q=vl.q;
                        int pp=ps*q; //p[ijk]
                        c.output_custom(format,id[ijk][q],p[ijk][pp],p[ijk][pp+1],p[ijk][pp+2],Config.default_radius,fp);
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
    void print_custom(String format,OutputStream fp) {
        CLoopAllPeriodic vl = new CLoopAllPeriodic(this);
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
     * vector.
     * \param[in] (x,y,z) the vector to test.
     * \param[out] (rx,ry,rz) the position of the particle whose Voronoi cell
     *                        contains the vector. This may point to a particle in
     *                        a periodic image of the primary domain.
     * \param[out] pid the ID of the particle.
     * \return True if a particle was found. If the container has no particles,
     * then the search will not find a Voronoi cell and false is returned. */
    public boolean find_voronoi_cell(double x,double y,double z,double[] rx,double[] ry,double[] rz,int[] pid) {
        int[] ai = new int[1],aj = new int[1],ak = new int[1],ci = new int[1],cj = new int[1],ck = new int[1],ijk = new int[1];
        ParticleRecord w = new ParticleRecord();
        double[] mrs = new double[1];

        // Remap the vector into the primary domain and then search for the
        // Voronoi cell that it is within
        double[] xout = new double[]{x}, yout = new double[]{y}, zout = new double[]{z};
        remap(ai,aj,ak,ci,cj,ck,xout,yout,zout,ijk);
        x = xout[0];
        y = yout[0];
        z = zout[0];
        vc.find_voronoi_cell(x,y,z,ci[0],cj[0],ck[0],ijk[0],w,mrs);

        if(w.ijk!=-1) {

            // Assemble the position vector of the particle to be returned,
            // applying a periodic remapping if necessary
            ci[0]+=w.di;
            if(ci[0]<0||ci[0]>=nx) ai[0]+=step_div(ci[0],nx);
            rx[0]=p[w.ijk][3*w.l]+ak[0]*bxz+aj[0]*bxy+ai[0]*bx;
            ry[0]=p[w.ijk][3*w.l+1]+ak[0]*byz+aj[0]*by;
            rz[0]=p[w.ijk][3*w.l+2]+ak[0]*bz;
            pid[0]=id[w.ijk][w.l];
            return true;
        }
        return false;

    }

    /** Computes the Voronoi cell for a particle currently being
     * referenced by a loop class.
     * \param[out] c a Voronoi cell class in which to store the
     * 		 computed cell.
     * \param[in] vl the loop class to use.
     * \return True if the cell was computed. If the cell cannot be
     * computed because it was removed entirely for some reason,
     * then the routine returns false. */
    public boolean compute_cell(VoronoiCellBase c,CLoopBase vl) {
        return vc.compute_cell(c,vl.ijk,vl.q,vl.i,vl.j,vl.k);
    }
    /** Computes the Voronoi cell for given particle.
     * \param[out] c a Voronoi cell class in which to store the
     * 		 computed cell.
     * \param[in] ijk the block that the particle is within.
     * \param[in] q the index of the particle within the block.
     * \return True if the cell was computed. If the cell cannot be
     * computed because it was removed entirely for some reason,
     * then the routine returns false. */
    public boolean compute_cell(VoronoiCellBase c,int ijk,int q) {
        int k = ijk/(nx*oy);
        int ijkt = ijk-(nx*oy)*k;
        int j = ijkt/nx;
        int i = ijkt-j*nx;
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
    public boolean compute_ghost_cell(VoronoiCellBase c,double x,double y,double z) {
        int[] ijk = new int[1];
        double[] xout = new double[]{x};
        double[] yout = new double[]{y};
        double[] zout = new double[]{z};
        put_locate_block(ijk,xout,yout,zout);
        x = xout[0];
        y = yout[0];
        z = zout[0];
        int pp=3*co[ijk[0]]; // p[ijk[0]]
        co[ijk[0]]++;
        p[ijk[0]][pp] = x;
        p[ijk[0]][pp+1] = y;
        p[ijk[0]][pp+2] = z;
        boolean q=compute_cell(c,ijk[0],co[ijk[0]]-1);
        co[ijk[0]]--;
        return q;
    }
    private VoroCompute vc;
}
