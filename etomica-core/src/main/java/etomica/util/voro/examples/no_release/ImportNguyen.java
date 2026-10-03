/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.util.voro.examples.no_release;

import etomica.util.collections.DoubleArrayList;
import etomica.util.voro.CLoopAll;
import etomica.util.voro.Container;
import etomica.util.voro.VoronoiCell;
import etomica.util.voro.examples.ResourceHelper;
import etomica.util.voro.examples.basic.Import;

import java.io.InputStream;

public class ImportNguyen {
    
    public static void main(String[] args) {
        // Set up finalants for the container geometry
        final double x_min=-5,x_max=5;
        final double y_min=-5,y_max=5;
        final double z_min=-5,z_max=5;

        // Set up the number of blocks that the container is divided into
        final int n_x=6,n_y=6,n_z=6;

        // Construct container
        Container con = new Container(-5,5,-5,5,0,10,6,6,6,false,false,false,8);

        // Import particles
        InputStream in = ResourceHelper.getStreamForFile("pack_ten_cube", Import.class);
        con.import_(in);

        // Loop over all the particles and compute the Voronoi cell for each
        int i;
        int id;
        double[] x = new double[1];
        double[] y = new double[1];
        double[] z = new double[1];
        DoubleArrayList vd = new DoubleArrayList();
        VoronoiCell c = new VoronoiCell();
        CLoopAll cl = new CLoopAll(con);
        if(cl.start()) do if(con.compute_cell(c,cl)) {

            // Get particle position and ID
            cl.pos(x,y,z);id=cl.pid();

            // Get face areas
            c.face_areas(vd);

            // Output information (additional diagnostics could be done
            // here)
            System.out.printf("ID %d (%.3f,%.3f,%.3f) :",id,x[0],y[0],z[0]);
            for(i=0;i<vd.size();i++) System.out.printf(" %.3f",vd.getDouble(i));
            System.out.println();
        } while (cl.inc());


    }
}
