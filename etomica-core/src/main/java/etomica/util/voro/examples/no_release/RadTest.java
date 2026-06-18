// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro.examples.no_release;

import etomica.util.voro.*;

public class RadTest {
    
    public static void main(String[] args) {
        // Set up constants for the container geometry
        final double x_min = -1, x_max = 1;
        final double y_min = -1, y_max = 1;
        final double z_min = -1, z_max = 1;
        final double cvol = (x_max - x_min) * (y_max - y_min) * (x_max - x_min);

        // Set up the number of blocks that the container is divided into
        final int n_x = 3, n_y = 3, n_z = 3;

        // Set the number of particles that are going to be randomly introduced
        final int particles = 1000;

        int i;
        double[] pos = new double[particles * 4];

        for (i = 0; i < particles; i++) {
            pos[4 * i + 0] = x_min + Math.random() * (x_max - x_min);
            pos[4 * i + 1] = y_min + Math.random() * (y_max - y_min);
            pos[4 * i + 2] = z_min + Math.random() * (z_max - z_min);
            pos[4 * i + 3] = Math.random();
        }

        for (i = 0; i <= 200; i++) {
            double mul = i * 0.01, vol = 0;

            ContainerPoly con = new ContainerPoly(x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z,
                    false, false, false, 8);

            for (int j = 0; j < particles; j++) {
                double x = pos[4 * j + 0];
                double y = pos[4 * j + 1];
                double z = pos[4 * j + 2];
                double r = pos[4 * j + 3] * mul;
                con.put(j, x, y, z, r);
            }

            String buf = String.format("rad_test_out/fr%d.pov", i);
//		FILE *pf=safe_fopen(buf,"w");
            int j = 0;
            CLoopAll cl = new CLoopAll(con);
            VoronoiCell c = new VoronoiCell();
            cl.start();
            do {
                if (con.compute_cell(c, cl)) {
                    vol += c.volume();
                    double[] xout = new double[1];
                    double[] yout = new double[1];
                    double[] zout = new double[1];
                    cl.pos(xout, yout, zout);
//				c.draw_pov(x,y,z,pf);
                    j++;
                }
            } while (cl.inc());

            System.out.printf("%g %d %g %g\n", mul, j, vol, vol - 8);
            if (Math.abs(vol-8) > 1e-13) Common.voro_fatal_error("Incorrect volume", Config.Voropp.INTERNAL_ERROR);

        }
    }
}
