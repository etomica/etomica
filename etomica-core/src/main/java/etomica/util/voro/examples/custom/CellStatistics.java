/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.util.voro.examples.custom;

import etomica.util.voro.VoronoiCell;

public class CellStatistics {

    public static void main(String[] args) {
        VoronoiCell v = new VoronoiCell();

        // Initialize the Voronoi cell to be a cube of side length 2, centered
        // on the origin
        v.init(-1,1,-1,1,-1,1);

        // Remove one edge of the cell with a single plane cut
        v.plane(1,1,0,2);

        // Output the Voronoi cell to a file in gnuplot format
        v.draw_gnuplot(0,0,0,"simple_cell.gnu");

        // Output vertex-based statistics
        System.out.printf("Total vertices      : %d\n",v.p);
        System.out.printf("Vertex positions    : ");v.output_vertices();System.out.println();
        System.out.printf("Vertex orders       : ");v.output_vertex_orders();System.out.println();
        System.out.printf("Max rad. sq. vertex : %g\n\n",0.25*v.max_radius_squared());

        // Output edge-based statistics
        System.out.printf("Total edges         : %d\n",v.number_of_edges());
        System.out.printf("Total edge distance : %g\n",v.total_edge_distance());
        System.out.printf("Face perimeters     : ");v.output_face_perimeters();System.out.println("\n");

        // Output face-based statistics
        System.out.printf("Total faces         : %d\n",v.number_of_faces());
        System.out.printf("Surface area        : %g\n",v.surface_area());
        System.out.printf("Face freq. table    : ");v.output_face_freq_table();System.out.println();
        System.out.printf("Face orders         : ");v.output_face_orders();System.out.println();
        System.out.printf("Face areas          : ");v.output_face_areas();System.out.println();
        System.out.printf("Face normals        : ");v.output_normals();System.out.println();
        System.out.printf("Face vertices       : ");v.output_face_vertices();System.out.println("\n");

        // Output volume-based statistics
        double[] xyz = v.centroid();
        System.out.printf("Volume              : %g\n"+
                "Centroid vector     : (%g,%g,%g)\n",v.volume(),xyz[0],xyz[1],xyz[2]);

    }


}
