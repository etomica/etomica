/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.util.voro.examples.custom;

import etomica.util.voro.Container;
import etomica.util.voro.ContainerPoly;
import etomica.util.voro.examples.ResourceHelper;

import java.io.InputStream;

public class Radical {

    public static void main(String[] args) {
        // Set up constants for the container geometry
        final double x_min=-3,x_max=3;
        final double y_min=-3,y_max=3;
        final double z_min=0,z_max=6;

        // Set up the number of blocks that the container is divided
        // into.
        final int n_x=3,n_y=3,n_z=3;

        // Create a container with the geometry given above, and make it
        // non-periodic in each of the three coordinates. Allocate space for
        // eight particles within each computational block. Import
        // the monodisperse test packing and output the Voronoi
        // tessellation in gnuplot and POV-Ray formats.
        Container con = new Container(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
                false,false,false,8);
        InputStream in = ResourceHelper.getStreamForFile("pack_six_cube", Radical.class);
        con.import_(in);
        con.draw_cells_gnuplot("pack_six_cube.gnu");
        con.draw_cells_pov("pack_six_cube_v.pov");
        con.draw_particles_pov("pack_six_cube_p.pov");

        // Create a container for polydisperse particles using the same
        // geometry as above. Import the polydisperse test packing and
        // output the Voronoi radical tessellation in gnuplot and POV-Ray
        // formats.
        ContainerPoly conp = new ContainerPoly(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
                false,false,false,8);
        in = ResourceHelper.getStreamForFile("pack_six_cube_poly", Radical.class);
        conp.import_(in);

        conp.draw_cells_gnuplot("pack_six_cube_poly.gnu");
        conp.draw_cells_pov("pack_six_cube_poly_v.pov");
        conp.draw_particles_pov("pack_six_cube_poly_p.pov");

    }
}
