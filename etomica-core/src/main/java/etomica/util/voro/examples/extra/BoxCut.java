/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.util.voro.examples.extra;

import etomica.util.voro.VoronoiCell;

public class BoxCut {

    public static void main(String[] args) {
        // Parameters controlling the center of the test box
        final double cx=1.5,cy=1.5,cz=1.5;

        double x,y,z;
        VoronoiCell v = new VoronoiCell();

        // Initialize the Voronoi cell to be a cube of side length 16, centered
        // on the origin
        v.init(-8,8,-8,8,-8,8);

        // Cut by a grid of points in a box of width one, centered on
        // (cx,cy,cz)
        for(x=cx-0.5;x<cx+0.55;x+=0.1) for(y=cy-0.5;y<cy+0.55;y+=0.1)
            for(z=cz-0.5;z<cz+0.55;z+=0.1) v.plane(x,y,z);

        // Output the Voronoi cell in gnuplot format
        v.draw_gnuplot(0,0,0,"box_cut.gnu");

        // Now make a small file that contains the test box
        v.init(cx-0.5,cx+0.5,cy-0.5,cy+0.5,cz-0.5,cz+0.5);
        v.draw_gnuplot(0,0,0,"box_cut.points");


    }
}
