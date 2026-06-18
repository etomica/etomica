// Cylindrical wall example code
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro.examples.no_release;

import etomica.util.voro.*;

public class CylinderInv {

    public static class WallCylinderInv implements Wall {
        public WallCylinderInv(double xc_,double yc_,double zc_,double xa_,double ya_,double za_,double rc_) {
            this(xc_, yc_, zc_, xa_, ya_, za_, rc_, -99);
        }
        public WallCylinderInv(double xc_,double yc_,double zc_,double xa_,double ya_,double za_,double rc_,int w_id_) {
            w_id = w_id_;
            xc = xc_;
            yc = yc_;
            zc = zc_;
            xa = xa_;
            ya = ya_;
            za = za_;
            asi = 1/(xa_*xa_+ya_*ya_+za_*za_);
            rc = rc_;
        }
        public boolean point_inside(double x,double y,double z) {
            double xd=x-xc,yd=y-yc,zd=z-zc;
            double pa=(xd*xa+yd*ya+zd*za)*asi;
            xd-=xa*pa;yd-=ya*pa;zd-=za*pa;
            return xd*xd+yd*yd+zd*zd>rc*rc;
        }

        public boolean cut_cell_base(VoronoiCellBase c, double x, double y, double z) {
            double xd=x-xc,yd=y-yc,zd=z-zc;
            double pa=(xd*xa+yd*ya+zd*za)*asi;
            xd-=xa*pa;yd-=ya*pa;zd-=za*pa;
            pa=xd*xd+yd*yd+zd*zd;
            if(pa>1e-5) {
                pa=2*(Math.sqrt(pa)*rc-pa);
                return c.nplane(-xd,-yd,-zd,-pa,w_id);
            }
            return true;
        }

        public boolean cut_cell(VoronoiCell c, double x, double y, double z) {return cut_cell_base(c,x,y,z);}
        public boolean cut_cell(VoronoiCellNeighbor c, double x, double y, double z) {return cut_cell_base(c,x,y,z);}
		private final int w_id;
		private final double xc,yc,zc,xa,ya,za,asi,rc;
    };


    public static void main(String[] args) {
        // Set up constants for the container geometry
        final double x_min=-6,x_max=6;
        final double y_min=-6,y_max=6;
        final double z_min=-6,z_max=6;

        // Set the computational grid size
        final int n_x=6,n_y=6,n_z=6;

        int i=0;
        double x,y,z;

        // Create a container with the geometry given above, and make it
        // non-periodic in each of the three coordinates. Allocate space for
        // eight particles within each computational block.
        Container con = new Container(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
                false,false,false,8);

        // Add a cylindrical wall to the container
        WallCylinderInv cyl = new WallCylinderInv(0,0,0,0,1,0,4.2);
        con.add_wall(cyl);

        // Place particles in a regular grid within the frustum, for points
        // which are within the wall boundaries
        for(z=-5.5;z<6;z+=1) for(y=-5.5;y<6;y+=1) for(x=-5.5;x<6;x+=1)
            if (con.point_inside(x,y,z)) {con.put(i,x,y,z);i++;}

        // Output the particle positions in POV-Ray format
        con.draw_particles_pov("cylinder_inv_p.pov");

        // Output the Voronoi cells in POV-Ray format
        con.draw_cells_pov("cylinder_inv_v.pov");

        // Output the particle positions in gnuplot format
        con.draw_particles("cylinder_inv.par");

        // Output the Voronoi cells in gnuplot format
        con.draw_cells_gnuplot("cylinder_inv.gnu");


    }
}
