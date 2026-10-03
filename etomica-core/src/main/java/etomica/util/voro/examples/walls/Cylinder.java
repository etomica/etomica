// Cylindrical wall example code
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro.examples.walls;

import etomica.util.voro.Container;
import etomica.util.voro.WallCylinder;
import etomica.util.voro.examples.ResourceHelper;

import java.io.InputStream;

public class Cylinder {

    public static void main(String[] args) {
        // Set up constants for the container geometry
        final double x_min=-6.5,x_max=6.5;
        final double y_min=-6.5,y_max=6.5;
        final double z_min=0,z_max=18.5;

        // Set the computational grid size
        final int n_x=7,n_y=7,n_z=14;

        // Create a container with the geometry given above, and make it
        // non-periodic in each of the three coordinates. Allocate space for
        // eight particles within each computational block.
        Container con = new Container(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
                false,false,false,8);

        // Add a cylindrical wall to the container
        WallCylinder cyl = new WallCylinder(0,0,0,0,0,1,6);
        con.add_wall(cyl);

        // Import the particles from a file
        InputStream in = ResourceHelper.getStreamForFile("pack_cylinder", Cylinder.class);
        con.import_(in);

        // Output the particle positions in POV-Ray format
        con.draw_particles_pov("cylinder_p.pov");

        // Output the Voronoi cells in POV-Ray format
        con.draw_cells_pov("cylinder_v.pov");

    }
}
