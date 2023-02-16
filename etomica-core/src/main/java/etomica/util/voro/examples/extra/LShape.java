// Irregular packing example code
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro.examples.extra;

import etomica.util.voro.Container;
import etomica.util.voro.VoronoiCell;
import etomica.util.voro.VoronoiCellNeighbor;
import etomica.util.voro.Wall;

public class LShape {

    // Create a wall class that will initialize the Voronoi cell to fill the
// L-shaped domain
    public static class WallLShape implements Wall {
        public WallLShape() {
            v = new VoronoiCell();
            v.init_l_shape();
            v.draw_gnuplot(0,0,0,"l_shape_init.gnu");
        };

        public boolean point_inside(double x,double y,double z) {return true;}
        public boolean cut_cell(VoronoiCell c,double x,double y,double z) {

            // Set the cell to be equal to the L-shape
            c.equalOperator(v);
            c.translate(-x,-y,-z);

            // Set the tolerance to 100, to make the code search
            // for cases where non-convex cells are cut in multiple
            // places
            c.big_tol=100;
            return true;
        }
        public boolean cut_cell(VoronoiCellNeighbor c, double x, double y, double z) {

            // Set the cell to be equal to the L-shape
            c.equalOperator(v);
            c.translate(-x,-y,-z);

            // Set the tolerance to 100, to make the code search
            // for cases where non-convex cells are cut in multiple
            // places
            c.big_tol=100;
            return true;
        }

        private final VoronoiCell v;
    };


    public static void main(String[] args) {

        // Set the number of particles that are going to be randomly introduced
        final int particles=20;

        int i=0;
        double x,y,z;
        // Create a container
        Container con = new Container(-1,1,-1,1,-1,1,5,5,5,false,false,false,8);

        // Create the L-shape wall class and add it to the container
        WallLShape wls = new WallLShape();
        con.add_wall(wls);

        // Add particles, making sure not to place any outside of the L-shape
        while(i<particles) {
            x=2*Math.random()-1;
            y=2*Math.random()-1;
            if(x<0&&y>0) continue;
            z=2*Math.random()-1;
            con.put(i,x,y,z);
            i++;
        }

        // Check the Voronoi cell volume; it should be 6
        System.out.printf("Voronoi cell volume: %.8g\n",con.sum_cell_volumes());

        // Save the particles and Voronoi cells in gnuplot format
        con.draw_particles("l_shape_p.gnu");
        con.draw_cells_gnuplot("l_shape_v.gnu");

    }

}

