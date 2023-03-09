// Timing test example code
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro.examples.timing;

import etomica.util.voro.Container;

public class TimingTest {
    
    public static void main(String[] args) {

        // Set up constants for the container geometry
        final double x_min=-1,x_max=1;
        final double y_min=-1,y_max=1;
        final double z_min=-1,z_max=1;
        
        // Set up the number of blocks that the container is divided into. If the
        // preprocessor variable NNN hasn't been passed to the code, then initialize it
        // to a good value. Otherwise, use the value that has been passed.
        int NNN = 26;
        final int n_x=NNN,n_y=NNN,n_z=NNN;
        
        // Set the number of particles that are going to be randomly introduced
        final int particles=100000;

        int i;double x,y,z;

        // Create a container with the geometry given above, and make it
        // periodic in each of the three coordinates. Allocate space for eight
        // particles within each computational block.
        Container con = new Container(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
                true,true,true,8);

        //Randomly add particles into the container
        for(i=0;i<particles;i++) {
            x=x_min+Math.random()*(x_max-x_min);
            y=y_min+Math.random()*(y_max-y_min);
            z=z_min+Math.random()*(z_max-z_min);
            con.put(i,x,y,z);
        }

        // Store the initial clock time
        long start=System.nanoTime();

        // Carry out a dummy computation of all cells in the entire container
        con.compute_all_cells();

        // Calculate the elapsed time and print it
        long end=System.nanoTime();
        double runtime=(end-start)/1e9;
        System.out.printf("%g\n",runtime);

    }
}
