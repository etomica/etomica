// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro.examples.no_release;

import etomica.util.voro.CLoopAll;
import etomica.util.voro.Container;
import etomica.util.voro.VoronoiCell;

public class Minkowski {

    public static void main(String[] args) {
        // Set up constants for the container geometry
        final double x_min=-1,x_max=1;
        final double y_min=-1,y_max=1;
        final double z_min=-1,z_max=1;
        
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
            con.put(i,x,y,z);
        }

        double[] vo = new double[400];
        double[] ar = new double[400];
        double[] tvo = new double[1];
        double[] tar = new double[1];
        for(i=0;i<400;i++) vo[i]=0;
        for(i=0;i<400;i++) ar[i]=0;

        VoronoiCell c = new VoronoiCell();
        CLoopAll cl = new CLoopAll(con);
        if(cl.start()) do if(con.compute_cell(c,cl)) {
            for(i=0;i<400;i++) {
                c.minkowski(i*0.005,tvo,tar);
                vo[i]+=tvo[0];ar[i]+=tar[0];
            }
        } while(cl.inc());


        for(i=0;i<400;i++) System.out.printf("%g %g %g\n",i*0.005,vo[i],ar[i]);
    }
}
