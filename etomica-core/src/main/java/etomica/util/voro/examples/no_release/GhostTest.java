// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro.examples.no_release;

import etomica.util.voro.ContainerPeriodicPoly;
import etomica.util.voro.VoronoiCell;

import java.io.FileOutputStream;
import java.io.IOException;

public class GhostTest {
    
    public static void main(String[] args) {

        // Set up finalants for the container geometry
        final double x_min=-1,x_max=1;
        final double y_min=-1,y_max=1;
        final double z_min=-1,z_max=1;
        final double cvol=(x_max-x_min)*(y_max-y_min)*(x_max-x_min);

        // Set up the number of blocks that the container is divided into
        final int n_x=6,n_y=6,n_z=6;

        // Set the number of particles that are going to be randomly introduced
        final int particles=20;

        int i;
        double x,y,z;

        // Create a container with the geometry given above, and make it
        // non-periodic in each of the three coordinates. Allocate space for
        // eight particles within each computational block
        ContainerPeriodicPoly con = new ContainerPeriodicPoly(2,0.5,2,0,0,2,n_x,n_y,n_z,8);

        // Randomly add particles into the container
        for(i=0;i<4;i++) {
            x=x_min+Math.random()*(x_max-x_min);
            y=y_min+Math.random()*(y_max-y_min);
            z=z_min+Math.random()*(z_max-z_min);
            System.out.println("con.put("+i+","+x+","+y+",0,1);");
            con.put(i,x,y,0,1);
        }

        // Output the particle positions in gnuplot format
        con.draw_particles("ghost_test_p.gnu");

        // Output the Voronoi cells in gnuplot format
        con.draw_cells_gnuplot("ghost_test_v.gnu");

        try {
            // Open file for test ghost cell
            FileOutputStream fp = new FileOutputStream("ghost_test_c.gnu");
            VoronoiCell c = new VoronoiCell();

            // Compute a range of ghost cells
//	for(y=-3.5;y<3.5;y+=0.05) if(con.compute_ghost_cell(c,1,y,0,1))
//		c.draw_gnuplot(1,y,0,fp);

            // Compute a single ghost cell
            if (con.compute_ghost_cell(c, 1.56, 0.67, 0, 1)) c.draw_gnuplot(1.56, 0.67, 0, fp);

            // Close ghost cell file
            fp.close();

            // Draw the domain
            con.draw_domain_gnuplot("ghost_test_d.gnu");
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }

    }
}
