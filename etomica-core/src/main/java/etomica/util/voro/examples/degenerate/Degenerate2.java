/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.util.voro.examples.degenerate;

import etomica.util.voro.VoronoiCell;

public class Degenerate2 {
    
    public static void main(String[] args) {

        // The total number of points to create as degenerate vertices
        final int points=100;

        // The number of planes that will be cut around each point
        final int N=64;
        final double step=2*Math.PI/N;

        // The angle (in radians) of the cutting planes from horizontal
        final double theta=0.04;

        double x,y,z,rsq,r;
        VoronoiCell v = new VoronoiCell();

        // Initialize the Voronoi cell to be a cube of side length 2, centered on
        // the origin
        v.init(-1,1,-1,1,-1,1);
        
        // Plane cutting
        for (int n=0; n<points; n++) {

            // Choose a random point
            x=2*Math.random()-1;
            y=2*Math.random()-1;
            z=2*Math.random()-1;

            // Skip it if it's outside the unit sphere or too close to the
            // origin
            rsq=x*x+y*y+z*z;
            if(rsq<0.01||rsq>1) {n--; continue;}

            // Rescale the point so that it has modulus 1, and then apply
            // plane cuts around this point
            r=1/ Math.sqrt(rsq);x*=r;y*=r;z*=r;
            rsq= Math.sqrt(x*x+y*y);r=z/rsq;
            double phi1 = Math.random()*step;
            for(int i=0; phi1+i*step<2*Math.PI;i++) {
                double phi = phi1 + i*step;
                v.plane(x * Math.cos(theta) + Math.sin(theta) * (-y * Math.cos(phi) / rsq - x * r * Math.sin(phi)),
                        y * Math.cos(theta) + Math.sin(theta) * (x * Math.cos(phi) / rsq - y * r * Math.sin(phi)),
                        z * Math.cos(theta) + Math.sin(theta) * rsq * Math.sin(phi), 1);
            }
        }

        // Output the Voronoi cell to a file in Gnuplot format
        v.draw_gnuplot(0,0,0,"degenerate2.gnu");

        // Optional POV-Ray output
        v.draw_pov(0,0,0,"degenerate2_v.pov");
        v.draw_pov_mesh(0,0,0,"degenerate2_m.pov");

    }
}
