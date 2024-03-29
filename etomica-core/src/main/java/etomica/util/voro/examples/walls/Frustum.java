// Frustum example code
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro.examples.walls;

import etomica.util.voro.Container;
import etomica.util.voro.WallCone;

public class Frustum {

    public static void main(String[] args) {
        int i = 0;

        // Create a container with the geometry given above, and make it
        // non-periodic in each of the three coordinates. Allocate space for
        // eight particles within each computational block.
        Container con = new Container(-1.2, 1.2, -1.2, 1.2, 0, 1, 14, 14, 7,
                false, false, false, 8);

        // Add a cylindrical wall to the container
        WallCone cone = new WallCone(0, 0, 2, 0, 0, -1, Math.atan(0.5));
        con.add_wall(cone);

        // Place particles in a regular grid within the frustum, for points
        // which are within the wall boundaries
        for (double z = 0.1;z<1;z+=0.2) {
            for (double y = -0.85; y < 1; y += 0.2) {
                for (double x = -0.95; x < 1; x += 0.2) {
                    if (con.point_inside(x, y, z)) {
                        con.put(i, x, y, z);
                        i++;
                    }
                }
            }
        }

        // Output the particle positions and Voronoi cells in Gnuplot format
        con.draw_particles("frustum_p.gnu");
        con.draw_cells_gnuplot("frustum_v.gnu");

        // Output the particle positions and Voronoi cells in POV-Ray format
        con.draw_particles_pov("frustum_p.pov");
        con.draw_cells_pov("frustum_v.pov");

        // Compute the volume of the Voronoi cells and compare it to the
        // exact frustum volume
        double evol = Math.PI * 1 * (0.5 * 0.5 + 0.5 * 1 + 1 * 1) / 3;
        double vvol = con.sum_cell_volumes();

        System.out.printf("Exact frustum volume : %g\n"+
                       "Voronoi cell volume  : %g\n"+
                       "Difference           : %g\n",evol, vvol, vvol-evol);
    }

}
