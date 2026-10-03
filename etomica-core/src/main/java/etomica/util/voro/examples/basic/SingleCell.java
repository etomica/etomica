// Single Voronoi cell example code
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro.examples.basic;

import etomica.util.voro.VoronoiCell;

public class SingleCell {

    public static void main(String[] args) {
        VoronoiCell v = new VoronoiCell();

        // Initialize the Voronoi cell to be a cube of side length 2, centered
        // on the origin
        v.init(-1,1,-1,1,-1,1);

        // Cut the cell by 250 random planes which are all a distance 1 away
        // from the origin, to make an approximation to a sphere
        for(int i=0;i<250;i++) {
            double x=2*Math.random()-1;
            double y=2*Math.random()-1;
            double z=2*Math.random()-1;
            double rsq=x*x+y*y+z*z;
            if(rsq>0.01&&rsq<1) {
                double r=1/Math.sqrt(rsq);
                x*=r;y*=r;z*=r;
                v.plane(x,y,z,1);
            }
        }

        // Output the Voronoi cell to a file, in the gnuplot format
        v.draw_gnuplot(0,0,0,"single_cell.gnu");
    }
}
