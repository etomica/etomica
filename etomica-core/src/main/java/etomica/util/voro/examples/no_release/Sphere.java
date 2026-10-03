// Sphere example code
//
// Author   : Chris H. Rycroft (Harvard SEAS)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro.examples.no_release;

import etomica.util.voro.Container;
import etomica.util.voro.WallSphere;

public class Sphere {

    public static void main(String[] args) {
        int i=0;
        double x,y,z,evol,vvol;

        // Create a container with the geometry given above, and make it
        // non-periodic in each of the three coordinates. Allocate space for
        // eight particles within each computational block.
        Container con = new Container(-5,5,-5,5,-5,5,6,6,6,
                false,false,false,8);

        // Add a cylindrical wall to the container
        WallSphere sph = new WallSphere(0,0,0,4);
        con.add_wall(sph);

        // Place particles in a regular grid within the frustum, for points
        // which are within the wall boundaries
        for(z=-4.5;z<5;z+=1) for(y=-4.5;y<5;y+=1) for(x=-4.5;x<5;x+=1) {
            if (con.point_inside(x,y,z)) {
                con.put(i,x,y,z);i++;
            }
        }

        // Output the particle positions and Voronoi cells in Gnuplot format
        con.draw_particles("sphere_p.gnu");
        con.draw_cells_gnuplot("sphere_v.gnu");

        // Compute the volume of the Voronoi cells and compare it to the
        // exact frustum volume
        evol=4/3.0*Math.PI*4*4*4;
        vvol=con.sum_cell_volumes();
        System.out.printf("Exact sphere volume : %g\n"+
                "Voronoi cell volume  : %g\n"+
                "Difference           : %g\n",evol,vvol,vvol-evol);

    }
}
