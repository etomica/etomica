// Splitting a Voronoi cell example code
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro.examples.no_release;

import etomica.util.voro.VoronoiCell;

public class SplitCell {
    
    public static void main(String[] args) {
        double x,y,z,rsq,r;
        VoronoiCell v = new VoronoiCell();
        VoronoiCell v2 = new VoronoiCell();

        // Initialize the Voronoi cell to be a cube of side length 2, centered
        // on the origin
        v.init(-1,1,-1,1,-1,1);

        // Cut the cell by 250 random planes which are all a distance 1 away
        // from the origin, to make an approximation to a sphere
        for(int i=0;i<250;i++) {
            x=2*Math.random()-1;
            y=2*Math.random()-1;
            z=2*Math.random()-1;
            rsq=x*x+y*y+z*z;
            if(rsq>0.01&&rsq<1) {
                r=1/Math.sqrt(rsq);x*=r;y*=r;z*=r;
                v.plane(x,y,z,1);
            }
        }

        // Make copy of the Voronoi cell
        v2.equalOperator(v);

        // Cut one copy in one direction, and the other copy in the other direction
        v.plane(1,0,0,0);
        v2.plane(-1,0,0,0);

        // Output the Voronoi cell to a file, in the gnuplot format
        v.draw_gnuplot(-0.05,0,0,"split_cell1.gnu");
        v2.draw_gnuplot(0.05,0,0,"split_cell2.gnu");

    }
}
