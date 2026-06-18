// Single Voronoi cell example code
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro.examples.no_release;

import etomica.util.voro.VoronoiCell;

public class SingleCell2D {

    public static void main(String[] args) {
        double x,y,rsq,r;
        VoronoiCell v = new VoronoiCell();

        // Initialize the Voronoi cell to be a square of side length 2 in the xy-plane. Set the
        // cell to be 1 high in the z-direction.
        v.init(-1,1,-1,1,-0.5,0.5);

        // Cut the cell by 250 random planes which are all a distance 1 away
        // from the origin, to make an approximation to a sphere
        for(int i=0;i<25;i++) {
            x=2*Math.random()-1;
            y=2*Math.random()-1;
            rsq=x*x+y*y;
            if(rsq>0.01&&rsq<1) {
                r=1/Math.sqrt(rsq);x*=r;y*=r;
                v.plane(x,y,0,1);
            }
        }

        // Output the Voronoi cell to a file, in the gnuplot format
        v.draw_gnuplot(0,0,0,"single_cell_2d.gnu");

    }
}
