/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.util.voro.examples.no_release;

import etomica.util.collections.IntArrayList;
import etomica.util.voro.CLoopAll;
import etomica.util.voro.Container;
import etomica.util.voro.VoronoiCellNeighbor;

public class NeighTest {
    
    public static void main(String[] args) {

        // Set up constants for the container geometry
        final double x_min=-1,x_max=1;
        final double y_min=-1,y_max=1;
        final double z_min=-1,z_max=1;
        final double cvol=(x_max-x_min)*(y_max-y_min)*(x_max-x_min);
        
        // Set up the number of blocks that the container is divided into
        final int n_x=6,n_y=6,n_z=6;
        
        // Set the number of particles that are going to be randomly introduced
        final int particles=40;

        int i;
        double x,y,z;
        VoronoiCellNeighbor c = new VoronoiCellNeighbor();
        IntArrayList neigh = new IntArrayList();

        // Create a container with the geometry given above, and make it
        // non-periodic in each of the three coordinates. Allocate space for
        // eight particles within each computational block
        Container con = new Container(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
                true,true,true,8);

        // Randomly add particles into the container
        for(i=0;i<particles;i++) {
            x=x_min+Math.random()*(x_max-x_min);
            y=y_min+Math.random()*(y_max-y_min);
            z=z_min+Math.random()*(z_max-z_min);
            con.put(i,x,y,z);
        }

        // Output the particle positions in gnuplot format
        con.draw_particles("random_points_p.gnu");

        // Output the Voronoi cells in gnuplot format
        con.draw_cells_gnuplot("random_points_v.gnu");

        // Loop over all of the particles and compute their Voronoi cells
        CLoopAll cl = new CLoopAll(con);
        if(cl.start()) do if(con.compute_cell(c,cl)) {
            double[] xout = new double[1];
            double[] yout = new double[1];
            double[] zout = new double[1];
            cl.pos(xout,yout,zout);
            x = xout[0];
            y = yout[0];
            z = zout[0];
            i=cl.pid();
            c.neighbors(neigh);

            // Print out the information about neighbors
            System.out.printf("Particle %2d at (% 1.3f,% 1.3f,% 1.3f):",i,x,y,z);
            for(int j=0;j<neigh.size();j++) System.out.printf(" %d",neigh.getInt(j));
            System.out.println();
        } while (cl.inc());

    }
}
