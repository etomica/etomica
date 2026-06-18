// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro.examples.basic;

import etomica.util.voro.Container;

public class RandomPoints {

    public static void main(String[] args) {

        final double x_min=-1,x_max=1;
        final double y_min=-1,y_max=1;
        final double z_min=-1,z_max=1;
        final double cvol=(x_max-x_min)*(y_max-y_min)*(z_max-z_min);

        // Set up the number of blocks that the container is divided into
        final int n_x=6,n_y=6,n_z=6;

        // Set the number of particles that are going to be randomly introduced
        final int particles=20;

        int i;
        double x,y,z;

        // Create a container with the geometry given above, and make it
        // non-periodic in each of the three coordinates. Allocate space for
        // eight particles within each computational block
        Container con = new Container(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
                false,false,false,8);

        // Randomly add particles into the container
        for(i=0;i<particles;i++) {
            x=x_min+Math.random()*(x_max-x_min);
            y=y_min+Math.random()*(y_max-y_min);
            z=z_min+Math.random()*(z_max-z_min);
            System.out.println(x+", "+y+", "+z);
            con.put(i,x,y,z);
        }

        // Sum up the volumes, and check that this matches the container volume
        double vvol=con.sum_cell_volumes();
        System.out.printf("Container volume : %g\n"+
                "Voronoi volume   : %g\n"+
                "Difference       : %g\n",cvol,vvol,vvol-cvol);

        // Output the particle positions in gnuplot format
        con.draw_particles("random_points_p.gnu");

        // Output the Voronoi cells in gnuplot format
        con.draw_cells_gnuplot("random_points_v.gnu");

    }
}
