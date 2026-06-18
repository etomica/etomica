// Example code demonstrating find_voronoi_cell function
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro.examples.interface_;

import etomica.util.voro.CLoopAll;
import etomica.util.voro.Container;
import etomica.util.voro.VoronoiCell;

import java.io.FileOutputStream;
import java.io.IOException;

public class FindVoroCell {

    public static void main(String[] args) {
        // The sampling distance for the grids of find_voronoi_cell calls
        final double h=0.05;

        // The cube of the sampling distance, corresponding the amount of volume
        // associated with a sample point
        final double hcube=h*h*h;

        // Set the number of particles that are going to be randomly introduced
        final int particles=20;

        int i;
        double x=0,y=0,z=0;
        double[] r = new double[1];
        double[] rx = new double[1];
        double[] ry = new double[1];
        double[] rz = new double[1];

        // Create a container with the geometry given above, and make it
        // non-periodic in each of the three coordinates. Allocate space for
        // eight particles within each computational block
        Container con = new Container(0,1,0,1,0,1,5,5,5,false,false,false,8);

        // Randomly add particles into the container
        for(i=0;i<particles;i++) {
            x=Math.random();
            y=Math.random();
            z=Math.random();
            con.put(i,x,y,z);
        }

        // Output the particle positions in gnuplot format
        con.draw_particles("find_voro_cell_p.gnu");

        // Scan a 2D slice in the container, and for each point in the slice,
        // find the Voronoi cell that the point is in. Store a vector
        int[] pid = new int[]{i};
        try {
            FileOutputStream f1 = new FileOutputStream("find_voro_cell.vec");
            for(x=0.5*h;x<1;x+=h) for(y=0.5*h;y<1;y+=h) {
                if(con.find_voronoi_cell(x,y,0.5,rx,ry,rz,pid))
                    f1.write(String.format("%g %g %g %g %g %g %g\n",x,y,0.5,rx[0]-x,ry[0]-y,rz[0]-0.5,
                            Math.sqrt((rx[0]-x)*(rx[0]-x)+(ry[0]-y)*(ry[0]-y)+(rz[0]-0.5)*(rz[0]-0.5))).getBytes());
                else {
                    System.err.printf("# find_voronoi_cell error for %g %g 0.5\n",x,y);
                }
            }
            f1.close();
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }


        // Create a blank array for storing the sampled Voronoi volumes
        int[] samp_v = new int[particles];
        for(i=0;i<particles;i++) samp_v[i]=0;

        // Scan over a grid covering the entire container, finding which
        // Voronoi cell each point is in, and tallying the result as a method
        // of sampling the volume of each Voronoi cell
        for(z=0.5*h;z<1;z+=h) for(y=0.5*h;y<1;y+=h) for(x=0.5*h;x<1;x+=h) {
            if(con.find_voronoi_cell(x,y,z,rx,ry,rz,pid)) samp_v[pid[0]]++;
            else System.err.printf("# find_voronoi_cell error for %g %g %g\n",x,y,z);
        }

        // Output the Voronoi cells in gnuplot format and a file with the
        // comparisons between the Voronoi cell volumes and the sampled volumes
        try {
            FileOutputStream f1 = new FileOutputStream("find_voro_cell.vol");
            FileOutputStream f2 = new FileOutputStream("find_voro_cell_v.gnu");
            CLoopAll cla = new CLoopAll(con);
            VoronoiCell c = new VoronoiCell();
            if (cla.start()) do if (con.compute_cell(c, cla)) {

                // Get the position and ID information for the particle
                // currently being considered by the loop. Ignore the radius
                // information.
                double[] xout = new double[]{x};
                double[] yout = new double[]{y};
                double[] zout = new double[]{z};
                cla.pos(pid, xout, yout, zout, r);
                x = xout[0];
                y = yout[0];
                z = zout[0];

                // Save and entry to the .vol file, storing both the computed
                // Voronoi cell volume, and the sampled volume based on the
                // number of grid points that were inside the cell
                f1.write(String.format("%d %g %g %g %g %g\n", pid[0], x, y, z, c.volume(), samp_v[pid[0]] * hcube).getBytes());

                // Draw the Voronoi cell
                c.draw_gnuplot(x, y, z, f2);
            } while (cla.inc());
            f1.close();
            f2.close();
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }
}
