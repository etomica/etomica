// Degenerate Voronoi cell example code
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro.examples.degenerate;

import etomica.util.voro.VoronoiCell;

public class Degenerate {

    public static void main(String[] args) {

        // The number of planes to be cut around each coordinate axis
        final int n=32;
        final double step=2*Math.PI/n;

        // The angle (in radians) of the cutting planes from horizontal
        final double theta=Math.PI/4-0.25;

        double x,y,z,phi;
        VoronoiCell v = new VoronoiCell();

        // Initialize the Voronoi cell to be a cube of side length 2, centered
        // on the origin
        v.init(-1,1,-1,1,-1,1);

        // Plane cutting
        for(phi=0;phi<2*Math.PI-0.5*step;phi+=step) {
            x=Math.cos(theta);y=Math.cos(phi)*Math.sin(theta);z=Math.sin(phi)*Math.sin(theta);
            System.out.println("phi "+phi+" "+x+" "+y+" "+z);
            v.plane(x,y,z,1);
            v.plane(-x,y,z,1);
            v.plane(y,x,z,1);
            v.plane(y,-x,z,1);
            v.plane(y,z,x,1);
            v.plane(y,z,-x,1);
        }

        // Check that the relation table is correct, and that there are no
        // duplicate edges
        v.check_relations();
        v.check_duplicates();

        // Output the Voronoi cell to a file in Gnuplot format
        v.draw_gnuplot(0,0,0,"degenerate.gnu");

        // Output the Voronoi cell to a file in POV-Ray format
        v.draw_pov(0,0,0,"degenerate_v.pov");

    }

}
