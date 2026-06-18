// Voronoi method to generate nanocrystalline grain boundaries
// Oct 18, 2011

package etomica.util.voro.examples.no_release;

import etomica.util.voro.CLoopOrder;
import etomica.util.voro.Container;
import etomica.util.voro.ParticleOrder;

import java.io.FileOutputStream;
import java.io.IOException;

public class ImportRahman {
    
    public static void main(String[] args) {

        // Box geometry
        final double x_min=-10,x_max=10;
        final double y_min=-10,y_max=10;
        final double z_min=-10,z_max=10;
        final double cvol=(x_max-x_min)*(y_max-y_min)*(x_max-x_min);

        // Number of blocks that the Box is divided into
        final int n_x=5,n_y=5,n_z=4;

        // Total no of particles

        final int particles=10000;

        int i;
        double x,y,z;

        // Creating Box and allcating 100 particles within each block

        Container con = new Container(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
                false,false,false,100);

        // Set up particle order class
        ParticleOrder po = new ParticleOrder();

        // Add particles into the Box	
        for(i=1;i<particles;i++) {
            x=x_min+Math.random()*(x_max-x_min);
            y=y_min+Math.random()*(y_max-y_min);
            z=z_min+Math.random()*(z_max-z_min);
            con.put(po,i,x,y,z);
        }

        // Setup an ordered loop class
        CLoopOrder cl = new CLoopOrder(con,po);

        try {
            // Customize output for LAMMPS, preserving ordering
            FileOutputStream fp = new FileOutputStream("lammps_input");
            con.print_custom(cl, "%i 1 %x %y %z", fp);
            fp.close();
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }

    }
}
