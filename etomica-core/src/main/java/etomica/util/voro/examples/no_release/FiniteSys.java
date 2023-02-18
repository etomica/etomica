// Irregular packing example code
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro.examples.no_release;

import etomica.util.collections.IntArrayList;
import etomica.util.voro.*;
import etomica.util.voro.examples.ResourceHelper;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;


public class FiniteSys {

    // ID for dodecahedron faces
    public static final int wid=-10;

    // Create a wall class that, whenever called, will replace the Voronoi cell
    // with a prescribed shape, in this case a dodecahedron
    public static class WallInitialShape implements Wall {


        public WallInitialShape() {
            final double Phi=0.5*(1+Math.sqrt(5.0));


            v = new VoronoiCellNeighbor();
            // Create a dodecahedron with neighbor information all
            // set to -10
            v.init(-2,2,-2,2,-2,2);
            v.nplane(0,Phi,1,wid);v.nplane(0,-Phi,1,wid);v.nplane(0,Phi,-1,wid);
            v.nplane(0,-Phi,-1,wid);v.nplane(1,0,Phi,wid);v.nplane(-1,0,Phi,wid);
            v.nplane(1,0,-Phi,wid);v.nplane(-1,0,-Phi,wid);v.nplane(Phi,1,0,wid);
            v.nplane(-Phi,1,0,wid);v.nplane(Phi,-1,0,wid);v.nplane(-Phi,-1,0,wid);
        };
        public boolean point_inside(double x,double y,double z) {return true;}
        public boolean cut_cell(VoronoiCell c, double x, double y, double z) {

            // Just ignore this case
            return true;
        }
        public boolean cut_cell(VoronoiCellNeighbor c,double x,double y,double z) {

            // Set the cell to be equal to the dodecahedron
            c.equalOperator(v);
            return true;
        }

        private final VoronoiCellNeighbor v;
    }

    // Determines whether any of the sides in the neighbor information are from the
    // initial dodecahedron
    public static boolean has_dodec_sides(IntArrayList vi) {
        for(int i=0;i<vi.size();i++) if(vi.getInt(i)==wid) return true;
        return false;
    }


    public static void main(String[] args) {
        // Set up finalants for the container geometry
        final double x_min=-15,x_max=15;
        final double y_min=-7,y_max=7;
        final double z_min=-15,z_max=15;

        // Golden ratio finalants
        final double phi=0.5*(1-Math.sqrt(5.0));

        // Set up the number of blocks that the container is divided
        // into.
        final int n_x=8,n_y=8,n_z=8;

        // Create a container with the geometry given above. This is bigger
        // than the particle packing itself.
        Container con = new Container(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
                false,false,false,8);

        // Create the "initial shape" wall class and add it to the container
        WallInitialShape wis = new WallInitialShape();
        con.add_wall(wis);

        // Import the irregular particle packing
        InputStream in = ResourceHelper.getStreamForFile("pack_semicircle", FiniteSys.class);
        con.import_(in);

        try {
            // Open files to save the "inside" and "outside" particles
            FileOutputStream finside = new FileOutputStream("finite_sys_in.pov");
            FileOutputStream foutside = new FileOutputStream("finite_sys_out.pov");

            // Loop over all particles
            double[] x = new double[1];
            double[] y = new double[1];
            double[] z = new double[1];
            IntArrayList vi = new IntArrayList();
            VoronoiCellNeighbor c = new VoronoiCellNeighbor();
            CLoopAll cl = new CLoopAll(con);
            if (cl.start()) do {

                // Get particle position
                cl.pos(x, y, z);

                // Remove half of the particles to see a cross-section
                if (y[0] < 0) continue;

                if (con.compute_cell(c, cl)) {

                    // Get the neighboring IDs of all the faces
                    c.neighbors(vi);

                    // Depending on whether any of the faces are
                    // from the original dodecahedron, print to
                    // the "inside" or "outside" file
                    FileOutputStream fp = has_dodec_sides(vi) ? foutside : finside;
                    fp.write(String.format("sphere{<%g,%g,%g>,s}\n", x[0], y[0], z[0]).getBytes());

                }
            } while (cl.inc());

            // Close the output files
            finside.close();
            foutside.close();
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }

    }
    
}
