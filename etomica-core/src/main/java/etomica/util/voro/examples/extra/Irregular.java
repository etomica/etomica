// Irregular packing example code
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro.examples.extra;

import etomica.util.voro.Container;
import etomica.util.voro.VoronoiCell;
import etomica.util.voro.VoronoiCellNeighbor;
import etomica.util.voro.Wall;
import etomica.util.voro.examples.ResourceHelper;

import java.io.InputStream;

public class Irregular {

    // Create a wall class that, whenever called, will replace the Voronoi cell
    // with a prescribed shape, in this case a dodecahedron
    public static class WallInitialShape implements Wall {

        // Golden ratio finalants
        static final double Phi=0.5*(1+Math.sqrt(5.0));

        public WallInitialShape() {
            v = new VoronoiCell();
            // Create a dodecahedron
            v.init(-2,2,-2,2,-2,2);
            v.plane(0,Phi,1);v.plane(0,-Phi,1);v.plane(0,Phi,-1);
            v.plane(0,-Phi,-1);v.plane(1,0,Phi);v.plane(-1,0,Phi);
            v.plane(1,0,-Phi);v.plane(-1,0,-Phi);v.plane(Phi,1,0);
            v.plane(-Phi,1,0);v.plane(Phi,-1,0);v.plane(-Phi,-1,0);
        };
        public boolean point_inside(double x,double y,double z) {return true;}
        public boolean cut_cell(VoronoiCell c,double x,double y,double z) {

            // Set the cell to be equal to the dodecahedron
            c.equalOperator(v);
            return true;
        }
        public boolean cut_cell(VoronoiCellNeighbor c, double x, double y, double z) {

            // Set the cell to be equal to the dodecahedron
            c.equalOperator(v);
            return true;
        }

        private final VoronoiCell v;

    }

    public static void main(String[] args) {

        // Set up finalants for the container geometry
        final double x_min=-6,x_max=6;
        final double y_min=-6,y_max=6;
        final double z_min=-3,z_max=9;

        // Set up the number of blocks that the container is divided
        // into.
        final int n_x=5,n_y=5,n_z=5;

        // Create a container with the geometry given above. This is bigger
        // than the particle packing itself.
        Container con = new Container(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
                false,false,false,8);

        // Create the "initial shape" wall class and add it to the container
        WallInitialShape wis = new WallInitialShape();
        con.add_wall(wis);

        // Import the irregular particle packing
        InputStream in = ResourceHelper.getStreamForFile("pack_irregular", Irregular.class);
        con.import_(in);

        // Save the particles and Voronoi cells in POV-Ray format
        con.draw_particles_pov("irregular_p.pov");
        con.draw_cells_pov("irregular_v.pov");
    }
}
