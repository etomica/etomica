// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro;

/** \brief A class representing a conical wall object.
 *
 * This class represents a cone wall object. */
public class WallCone implements Wall {

    public WallCone(double xc_,double yc_,double zc_,double xa_,double ya_,double za_,double ang) {
        this(xc_, yc_, zc_, xa_, ya_, za_, ang, -99);
    }

        /** Constructs a cone wall object.
         * \param[in] (xc_,yc_,zc_) the apex of the cone.
         * \param[in] (xa_,ya_,za_) a vector pointing along the axis of
         *			    the cone.
         * \param[in] ang the angle (in radians) of the cone, measured
         *		  from the axis.
         * \param[in] w_id_ an ID number to associate with the wall for
         *		    neighbor tracking. */
    public WallCone(double xc_,double yc_,double zc_,double xa_,double ya_,double za_,double ang,int w_id_) {
        w_id = w_id_;
        xc = xc_;
        yc = yc_;
        zc = zc_;
        xa = xa_;
        ya = ya_;
        za = za_;
        asi = 1/(xa_*xa_+ya_*ya_+za_*za_);
        gra = Math.tan(ang);
        sang = Math.sin(ang);
        cang = Math.cos(ang);
    }

    /** Tests to see whether a point is inside the cone wall object.
     * \param[in] (x,y,z) the vector to test.
     * \return True if the point is inside, false if the point is outside. */
    public boolean point_inside(double x,double y,double z) {
        double xd=x-xc,yd=y-yc,zd=z-zc,pa=(xd*xa+yd*ya+zd*za)*asi;
        xd-=xa*pa;yd-=ya*pa;zd-=za*pa;
        pa*=gra;
        if (pa<0) return false;
        pa*=pa;
        return xd*xd+yd*yd+zd*zd<pa;
    }

    /** Cuts a cell by the cone wall object. The conical wall is
     * approximated by a single plane applied at the point on the cone which is
     * closest to the center of the cell. This works well for particle arrangements
     * that are packed against the wall, but loses accuracy for sparse particle
     * distributions.
     * \param[in,out] c the Voronoi cell to be cut.
     * \param[in] (x,y,z) the location of the Voronoi cell.
     * \return True if the cell still exists, false if the cell is deleted. */
    public boolean cut_cell_base(VoronoiCellBase c,double x,double y,double z) {
        double xd=x-xc,yd=y-yc,zd=z-zc,xf,yf,zf,q,pa=(xd*xa+yd*ya+zd*za)*asi;
        xd-=xa*pa;yd-=ya*pa;zd-=za*pa;
        pa=xd*xd+yd*yd+zd*zd;
        if(pa>1e-5) {
            pa=1/Math.sqrt(pa);
            q=Math.sqrt(asi);
            xf=-sang*q*xa+cang*pa*xd;
            yf=-sang*q*ya+cang*pa*yd;
            zf=-sang*q*za+cang*pa*zd;
            pa=2*(xf*(xc-x)+yf*(yc-y)+zf*(zc-z));
            return c.nplane(xf,yf,zf,pa,w_id);
        }
        return true;
    }


    public boolean cut_cell(VoronoiCell c,double x,double y,double z) {return cut_cell_base(c,x,y,z);}
    public boolean cut_cell(VoronoiCellNeighbor c,double x,double y,double z) {return cut_cell_base(c,x,y,z);}

    private final int w_id;
    private final double xc,yc,zc,xa,ya,za,asi,gra,sang,cang;

}
