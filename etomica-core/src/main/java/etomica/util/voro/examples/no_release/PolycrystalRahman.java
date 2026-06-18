/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.util.voro.examples.no_release;

import etomica.util.voro.CLoopAll;
import etomica.util.voro.Container;
import etomica.util.voro.VoronoiCell;

import java.io.FileOutputStream;
import java.io.IOException;

public class PolycrystalRahman {
    
    public static void main(String[] args) {

        // Box geometry
        final double x_min=0,x_max=81;
        final double y_min=0,y_max=81;
        final double z_min=0,z_max=81;
        final double cvol=(x_max-x_min)*(y_max-y_min)*(x_max-x_min);
        
        // Total number of particles
        final int particles=20;
        
        // Lattice size
        final double h=4;

        int j=0,c_id;
        double x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3,xt,yt,zt;
        double theta,cth,sth,r,v,vx,vy,vz;

        // Create the container class
        Container con = new Container(x_min,x_max,y_min,y_max,z_min,z_max,10,10,10,
                false,false,false,8);

        // Add generators to the box
        for(int i=1;i<=particles;i++) {
            x=x_min+Math.random()*(x_max-x_min);
            y=y_min+Math.random()*(y_max-y_min);
            z=z_min+Math.random()*(z_max-z_min);
            con.put(i,x,y,z);
        }

        // Open the file for the particle positions
        try {
            FileOutputStream fp = new FileOutputStream("lammps_input");

            // Create a loop class to iterate over all of the generators
            // in the box
            CLoopAll cl = new CLoopAll(con);
            VoronoiCell c = new VoronoiCell();
            if (cl.start()) do if (con.compute_cell(c, cl)) {

                // Generate the first vector of an orthonormal basis
                x1 = 2 * Math.random() - 1;
                y1 = 2 * Math.random() - 1;
                z1 = 2 * Math.random() - 1;
                r = 1 / Math.sqrt(x1 * x1 + y1 * y1 + z1 * z1);
                x1 *= r;
                y1 *= r;
                z1 *= r;

                // Construct a second perpendicular vector
                if (Math.abs(x1) > 0.5 || Math.abs(y1) > 0.5) {
                    r = 1 / Math.sqrt(x1 * x1 + y1 * y1);
                    x2 = -y1 * r;
                    y2 = x1 * r;
                    z2 = 0;
                } else {
                    r = 1 / Math.sqrt(x1 * x1 + z1 * z1);
                    x2 = -z1 * r;
                    y2 = 0;
                    z2 = x1 * r;
                }

                // Construct a third perpendicular vector using the vector product
                x = y2 * z1 - z2 * y1;
                y = z2 * x1 - x2 * z1;
                z = x2 * y1 - y2 * x1;

                // Take a random rotation of the second and third vectors
                theta = 2 * Math.PI * Math.random();

                cth = Math.cos(theta);
                sth = Math.sin(theta);
                x3 = x * cth + x2 * sth;
                x2 = -x * sth + x2 * cth;
                y3 = y * cth + y2 * sth;
                y2 = -y * sth + y2 * cth;
                z3 = z * cth + z2 * sth;
                z2 = -z * sth + z2 * cth;

                // Get a bound on how far to search
                r = Math.sqrt(0.25 * c.max_radius_squared());
                v = (int) ((r / h) + 2) * h;

                // Add small random displacement to lattice positioning,
                // so that it's not always perfectly aligned with the generator
                vx = -v + h * Math.random();
                vy = -v + h * Math.random();
                vz = -v + h * Math.random();

                // Print diagnostic information about this generator
                c_id = cl.pid();
                double[] xc = new double[1];
                double[] yc = new double[1];
                double[] zc = new double[1];
                cl.pos(xc, yc, zc);
                System.out.printf("Generator %d at (%g,%g,%g), random basis:\n", c_id, xc[0], yc[0], zc[0]);
                System.out.printf("%g %g %g\n", x1, y1, z1);
                System.out.printf("%g %g %g\n", x2, y2, z2);
                System.out.printf("%g %g %g\n\n", x3, y3, z3);

                // Loop over a local region of points
                for (z = vx; z <= v; z += h)
                    for (y = vy; y <= v; y += h)
                        for (x = vz; x <= v; x += h) {

                            // Construct points rotated into the random basis
                            xt = xc[0] + x * x1 + y * x2 + z * x3;
                            yt = yc[0] + x * y1 + y * y2 + z * y3;
                            zt = zc[0] + x * z1 + y * z2 + z * z3;

                            // Skip if this lies outside the container
                            if (xt < x_min || xt > x_max || yt < y_min || yt > y_max || zt < z_min || zt > z_max)
                                continue;

                            // Find the nearest generator
                            double[] rx = new double[1];
                            double[] ry = new double[1];
                            double[] rz = new double[1];
                            int[] i = new int[1];
                            con.find_voronoi_cell(xt, yt, zt, rx, ry, rz, i);

                            // If the nearest generator matches, then save this point
                            if (i[0] == c_id) {
                                fp.write(String.format("%d %g %g %g\n", j, xt, yt, zt).getBytes());
                                j++;
                            }
                        }
            } while (cl.inc());

            // Close the output file
            fp.close();
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }

        // Output files for diagnosic purposes
        con.draw_particles("lammps_generators");
        con.draw_cells_gnuplot("lammps_cells");

    }
}
