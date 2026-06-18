// Example code demonstrating the loop classes
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro.examples.interface_;

import etomica.util.voro.*;

import java.io.FileOutputStream;
import java.io.IOException;

public class Loops {
    
    public static void main(String[] args) {
        // Constants determining the configuration of the tori
        final double dis=1.25,mjrad=2.5,mirad=0.95,trad=mjrad+mirad;

        // Set the number of particles that are going to be randomly introduced
        final int particles=100000;

        int i;
        double r;
        VoronoiCell c = new VoronoiCell();

        // Create a container as a non-periodic 10 by 10 by 10 box
        Container con = new Container(-5,5,-5,5,-5,5,26,26,26,false,false,false,8);
        ParticleOrder po = new ParticleOrder();

        // Randomly add particles into the container
        for(i=0;i<particles;i++) {
            double x=10*Math.random()-5;
            double y=10*Math.random()-5;
            double z=10*Math.random()-5;

            // If the particle lies within the first torus, store it in the
            // ordering class when adding to the container
            r=Math.sqrt((x-dis)*(x-dis)+y*y);
            if((r-mjrad)*(r-mjrad)+z*z<mirad) con.put(po,i,x,y,z);
            else con.put(i,x,y,z);
        }

        // Compute Voronoi cells for the first torus. Here, the points
        // previously stored in the ordering class are looped over.
        try {
            FileOutputStream f1 = new FileOutputStream("loops1_m.pov");
            FileOutputStream f2 = new FileOutputStream("loops1_v.pov");
            CLoopOrder clo = new CLoopOrder(con, po);
            if (clo.start()) do if (con.compute_cell(c, clo)) {

                // Get the position of the current particle under consideration
                double[] x = new double[1];
                double[] y = new double[2];
                double[] z = new double[3];
                clo.pos(x, y, z);

                // Save a POV-Ray mesh to one file and a cylinder/sphere
                // representation to the other file
                c.draw_pov_mesh(x[0], y[0], z[0], f1);
                c.draw_pov(x[0], y[0], z[0], f2);
            } while (clo.inc());
            f1.close();
            f2.close();
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }

        try {
            // Compute Voronoi cells for the second torus. Here, the subset loop is
            // used to search over the blocks overlapping the torus, and then each
            // particle is individually tested.
            FileOutputStream f1 = new FileOutputStream("loops2_m.pov");
            FileOutputStream f2 = new FileOutputStream("loops2_v.pov");
            CLoopSubset cls = new CLoopSubset(con);
            cls.setup_box(-dis - trad, -dis + trad, -mirad, mirad, -trad, trad, false);
            if (cls.start()) do {

                // Get the position of the current particle under consideration
                double[] x = new double[1];
                double[] y = new double[1];
                double[] z = new double[1];
                cls.pos(x, y, z);

                // Test whether this point is within the torus, and if so,
                // compute and save the Voronoi cell
                r = Math.sqrt((x[0] + dis) * (x[0] + dis) + z[0] * z[0]);
                if ((r - mjrad) * (r - mjrad) + y[0] * y[0] < mirad && con.compute_cell(c, cls)) {
                    c.draw_pov_mesh(x[0], y[0], z[0], f1);
                    c.draw_pov(x[0], y[0], z[0], f2);
                }
            } while (cls.inc());
            f1.close();
            f2.close();
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }
}
