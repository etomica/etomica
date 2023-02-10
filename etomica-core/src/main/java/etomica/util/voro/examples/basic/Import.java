// File import example code
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro.examples.basic;

import etomica.util.voro.Container;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStream;

public class Import {

    public static void main(String[] args) {
        // Set up constants for the container geometry
        final double x_min=-5,x_max=5;
        final double y_min=-5,y_max=5;
        final double z_min=0,z_max=10;

        // Set up the number of blocks that the container is divided into
        final int n_x=6,n_y=6,n_z=6;

        // Create a container with the geometry given above, and make it
        // non-periodic in each of the three coordinates. Allocate space for
        // eight particles within each computational block
        Container con = new Container(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
                false,false,false,8);

        //Randomly add particles into the container
        InputStream in = null;
        String fn = "pack_ten_cube";
        if (new File(fn).exists()) {
            // local file exists; use that
            try {
                in = new FileInputStream(fn);
            } catch (FileNotFoundException e) {
                throw new RuntimeException(e);
            }
        }
        else {
            // use file from resources
            in = Import.class.getResourceAsStream(fn);
        }
        con.import_(in);

        // Save the Voronoi network of all the particles to text files
        // in gnuplot and POV-Ray formats
        con.draw_cells_gnuplot("pack_ten_cube.gnu");
        con.draw_cells_pov("pack_ten_cube_v.pov");

        // Output the particles in POV-Ray format
        con.draw_particles_pov("pack_ten_cube_p.pov");

        // use ouput pov files in conjunction with import.pov
    }
}
