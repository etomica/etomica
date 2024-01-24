// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro;

import java.io.*;

/** \brief A class for storing an arbitrary number of particles without radius
 * information, prior to setting up a container geometry.
 *
 * The pre_container class is an extension of the pre_container_base class for
 * cases when no particle radius information is available. */
public class PreContainer extends PreContainerBase {

    /** The class constructor sets up the geometry of container,
     * initializing the minimum and maximum coordinates in each
     * direction.
     * \param[in] (ax_,bx_) the minimum and maximum x coordinates.
     * \param[in] (ay_,by_) the minimum and maximum y coordinates.
     * \param[in] (az_,bz_) the minimum and maximum z coordinates.
     * \param[in] (xperiodic_,yperiodic_,zperiodic_ ) flags setting whether the
     *                                                container is periodic in
     *                                                each coordinate direction. */
    public PreContainer(double ax_,double bx_,double ay_,double by_,double az_,double bz_,
                  boolean xperiodic_,boolean yperiodic_,boolean zperiodic_) {
        super(ax_, bx_, ay_, by_, az_, bz_, xperiodic_, yperiodic_, zperiodic_, 3);
    }

    /** Stores a particle ID and position, allocating a new memory chunk if
     * necessary. For coordinate directions in which the container is not periodic,
     * the routine checks to make sure that the particle is within the container
     * bounds. If the particle is out of bounds, it is not stored.
     * \param[in] n the numerical ID of the inserted particle.
     * \param[in] (x,y,z) the position vector of the inserted particle. */
    public void put(int n,double x,double y,double z) {
        if((xperiodic||(x>=ax&&x<=bx))&&(yperiodic||(y>=ay&&y<=by))&&(zperiodic||(z>=az&&z<=bz))) {
            if(ch_id==e_id) new_chunk();
            pre_id[end_id][ch_id] = n;
            ch_id++;
            pre_p[end_p][ch_p+0] = x;
            pre_p[end_p][ch_p+1] = y;
            pre_p[end_p][ch_p+2] = z;
            ch_p+=3;
        }
        else if (Config.VOROPP_REPORT_OUT_OF_BOUNDS) {
            System.err.printf("Out of bounds: (x,y,z)=(%g,%g,%g)\n", x, y, z);
        }
    }

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
    /** Imports particles from a file.
     * \param[in] filename the name of the file to read from. */
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
    public void setup(Container con) {
        for (int c_id=0, c_p=0; c_id<end_id; c_id++, c_p++) {
            int ide=Config.pre_container_chunk_size;
            for (int pp=0, idp=0; idp<ide; idp++, pp+=3) {
                int n = pre_id[c_id][idp];
                double x = pre_p[c_p][pp+0];
                double y = pre_p[c_p][pp+1];
                double z = pre_p[c_p][pp+2];
                con.put(n,x,y,z);
                idp++;
            }
        }
        for (int idp=0, pp=0; idp<ch_id; idp++, pp+=3) {
            int n = pre_id[end_id][idp];
            double x = pre_p[end_id][pp+0];
            double y = pre_p[end_id][pp+1];
            double z = pre_p[end_id][pp+2];
            con.put(n,x,y,z);
        }
    }

    public void setup(ParticleOrder vo,Container con) {
        for (int c_id=0, c_p=0; c_id<end_id; c_id++, c_p++) {
            int ide=Config.pre_container_chunk_size;
            for (int pp=0, idp=0; idp<ide; idp++, pp++) {
                int n = pre_id[c_id][idp];
                double x = pre_p[c_p][pp+0];
                double y = pre_p[c_p][pp+1];
                double z = pre_p[c_p][pp+2];
                con.put(vo,n,x,y,z);
                idp++;
            }
        }
        for (int idp=0, pp=0; idp<ch_id; idp++, pp+=3) {
            int n = pre_id[end_id][idp];
            double x = pre_p[end_id][pp+0];
            double y = pre_p[end_id][pp+1];
            double z = pre_p[end_id][pp+2];
            con.put(vo,n,x,y,z);
        }
    }

}
