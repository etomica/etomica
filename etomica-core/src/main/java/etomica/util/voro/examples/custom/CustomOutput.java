/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.util.voro.examples.custom;

import etomica.util.voro.Container;
import etomica.util.voro.examples.ResourceHelper;

import java.io.InputStream;

public class CustomOutput {

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
        // eight particles within each computational block.
        Container con = new Container(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
                false,false,false,8);

        // Import the monodisperse test packing and output the Voronoi
        // tessellation in gnuplot and POV-Ray formats.
        InputStream in = ResourceHelper.getStreamForFile("pack_six_cube", CustomOutput.class);
        con.import_(in);

        // Do a custom output routine to store the number of vertices, edges,
        // and faces of each Voronoi cell
        con.print_custom(
                "ID=%i, pos=(%x,%y,%z), vertices=%w, edges=%g, faces=%s",
                "packing.custom1");

        // Do a custom output routine to store a variety of face-based
        // statistics. Store the particle ID and position, the number of faces
        // the total face area, the order of each face, the areas of each face,
        // the vertices making up each face, and the neighboring particle (or
        // wall) corresponding to each face.
        con.print_custom("%i %q %s %F %a %f %t %l %n","packing.custom2");

        // Do a custom output routine that outputs the particle IDs and
        // positions, plus the volume and the centroid position relative to the
        // particle center
        con.print_custom("%i %q %v %c","packing.custom3");

        // Also create POV-Ray output of the Voronoi cells for use in the
        // rendering
        con.draw_cells_pov("pack_six_cube_v.pov");

    }
}
