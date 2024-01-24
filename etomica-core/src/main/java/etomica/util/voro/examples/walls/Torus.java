// Custom wall class example code
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro.examples.walls;

import etomica.util.voro.*;
import etomica.util.voro.examples.ResourceHelper;

import java.io.InputStream;


public class Torus {

    // This class creates a custom toroidal wall object that is centered on the
    // origin and is aligned with the xy plane. It is derived from the pure virtual
    // "wall" class. The "wall" class contains virtual functions for cutting the
    // Voronoi cell in response to a wall, and for telling whether a given point is
    // inside the wall or not. In this derived class, specific implementations of
    // these functions are given.
    public static class WallTorus implements Wall {

        public WallTorus(double imjr, double imnr) {
            this(imjr, imnr, -99);
        }

        // The wall constructor initializes constants for the major and
        // minor axes of the torus. It also initializes the wall ID
        // number that is used when the plane cuts are made. This is
        // only tracked with the voronoicell_neighbor class and is
        // ignored otherwise. It can be omitted, and then an arbitrary
        // value of -99 is used.
        WallTorus(double imjr,double imnr,int iw_id) {
            w_id = iw_id;
            mjr = imjr;
            mnr = imnr;
        }

        // This returns true if a given vector is inside the torus, and
        // false if it is outside. For the current example, this
        // routine is not needed, but in general it would be, for use
        // with the point_inside() routine in the container class.
        public boolean point_inside(double x,double y,double z) {
            double temp=Math.sqrt(x*x+y*y)-mjr;
            return temp*temp+z*z<mnr*mnr;
        }

        // This template takes a reference to a voronoicell or
        // voronoicell_neighbor object for a particle at a vector
        // (x,y,z), and makes a plane cut to to the object to account
        // for the toroidal wall
        public boolean cut_cell_base(VoronoiCellBase c, double x, double y, double z) {
            double orad=Math.sqrt(x*x+y*y);
            double odis=orad-mjr;
            double ot=odis*odis+z*z;

            // Unless the particle is within 1% of the minor
            // radius, then a plane cut is made
            if(ot>0.01*mnr*mnr) {
                ot=2*mnr/Math.sqrt(ot)-2;
                z*=ot;
                odis*=ot/orad;
                x*=odis;
                y*=odis;
                return c.nplane(x,y,z,w_id);
            }
            return true;
        }

        // These virtual functions are called during the cell
        // computation in the container class. They call instances of
        // the template given above.
        public boolean cut_cell(VoronoiCell c, double x,
                                double y, double z) {return cut_cell_base(c,x,y,z);}
        public boolean cut_cell(VoronoiCellNeighbor c, double x,
                                double y, double z) {return cut_cell_base(c,x,y,z);}

        // The ID number associated with the wall
        private final int w_id;
        // The major radius of the torus
        private final double mjr;
        // The minor radius of the torus
        private final double mnr;
    };


    public static void main(String[] args) {
        // Major and minor torus radii
        final double arad=9,brad=3.5;
        
        // The outer radius of the torus, that determines how big the container should
        // be
        final double crad=arad+brad;
        
        // Set up constants for the container geometry
        final double x_min=-crad-0.5,x_max=crad+0.5;
        final double y_min=-crad-0.5,y_max=crad+0.5;
        final double z_min=-brad-0.5,z_max=brad+0.5;
        
        // Set the computational grid size
        final int n_x=10,n_y=10,n_z=3;

        // Create a container with the geometry given above, and make it
        // non-periodic in each of the three coordinates. Allocate space for
        // eight particles within each computational block.
        Container con = new Container(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
                false,false,false,8);

        // Add the custom toroidal wall to the container
        WallTorus tor = new WallTorus(arad,brad);
        con.add_wall(tor);

        // Import the particles from a file
        InputStream in = ResourceHelper.getStreamForFile("pack_torus", Torus.class);
        con.import_(in);

        // Output the particle positions in POV-Ray format
        con.draw_particles_pov("torus_p.pov");

        // Output the Voronoi cells in POV-Ray format
        con.draw_cells_pov("torus_v.pov");

    }
}
