// Tetrahedron example code
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro.examples.walls;

import etomica.util.voro.Container;
import etomica.util.voro.WallPlane;

public class Tetrahedron {
    
    public static void main(String[] args) {
        // Set up constants for the container geometry
        final double x_min=-2,x_max=2;
        final double y_min=-2,y_max=2;
        final double z_min=-2,z_max=2;
        
        // Set up the number of blocks that the container is divided
        // into
        final int n_x=7,n_y=7,n_z=7;
        
        // Set the number of particles that are going to be randomly
        // introduced
        final int particles=64;

        int i=0;
        double x,y,z;

        // Create a container with the geometry given above, and make it
        // non-periodic in each of the three coordinates. Allocate space for 8
        // particles within each computational block.
        Container con = new Container(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
                false,false,false,8);

        // Add four plane walls to the container to make a tetrahedron
        WallPlane p1 = new WallPlane(1,1,1,1);
        con.add_wall(p1);
        WallPlane p2 = new WallPlane(-1,-1,1,1);
        con.add_wall(p2);
        WallPlane p3 = new WallPlane(1,-1,-1,1);
        con.add_wall(p3);
        WallPlane p4 = new WallPlane(-1,1,-1,1);
        con.add_wall(p4);

        // Randomly insert particles into the container, checking that they lie
        // inside the tetrahedron
        while(i<particles) {
            x=x_min+Math.random()*(x_max-x_min);
            y=y_min+Math.random()*(y_max-y_min);
            z=z_min+Math.random()*(z_max-z_min);
            if (con.point_inside(x,y,z)) {
                con.put(i,x,y,z);i++;
            }
        }

        // Output the particle positions and the Voronoi cells in Gnuplot and
        // POV-Ray formats
        con.draw_particles("tetrahedron_p.gnu");
        con.draw_cells_gnuplot("tetrahedron_v.gnu");
        con.draw_particles_pov("tetrahedron_p.pov");
        con.draw_cells_pov("tetrahedron_v.pov");

    }
}
