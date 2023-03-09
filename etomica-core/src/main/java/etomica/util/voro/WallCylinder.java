// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro;

public class WallCylinder implements Wall {

    public WallCylinder(double xc_, double yc_, double zc_, double xa_, double ya_, double za_, double rc_) {
        this(xc_, yc_, zc_, xa_, ya_, za_, rc_, -99);
    }
    /** Constructs a cylinder wall object.
     * \param[in] (xc_,yc_,zc_) a point on the axis of the
     *			    cylinder.
     * \param[in] (xa_,ya_,za_) a vector pointing along the
     *			    direction of the cylinder.
     * \param[in] rc_ the radius of the cylinder
     * \param[in] w_id_ an ID number to associate with the wall for
     *		    neighbor tracking. */
    public WallCylinder(double xc_,double yc_,double zc_,double xa_,double ya_,double za_,double rc_,int w_id_) {
        xc = xc_;
        yc = yc_;
        zc = zc_;
        xa = xa_;
        ya = ya_;
        za = za_;
        rc = rc_;
        w_id = w_id_;
        asi = 1/(xa_*xa_ + ya_*ya_ + za_*za_);
    }

    /** Tests to see whether a point is inside the cylindrical wall object.
     * \param[in] (x,y,z) the vector to test.
     * \return True if the point is inside, false if the point is outside. */
    public boolean point_inside(double x,double y,double z) {
        double xd=x-xc,yd=y-yc,zd=z-zc;
        double pa=(xd*xa+yd*ya+zd*za)*asi;
        xd-=xa*pa;yd-=ya*pa;zd-=za*pa;
        return xd*xd+yd*yd+zd*zd<rc*rc;
    }

    /** Cuts a cell by the cylindrical wall object. The cylindrical wall is
     * approximated by a single plane applied at the point on the cylinder which is
     * closest to the center of the cell. This works well for particle arrangements
     * that are packed against the wall, but loses accuracy for sparse particle
     * distributions.
     * \param[in,out] c the Voronoi cell to be cut.
     * \param[in] (x,y,z) the location of the Voronoi cell.
     * \return True if the cell still exists, false if the cell is deleted. */
    public boolean cut_cell_base(VoronoiCellBase c,double x,double y,double z) {
        double xd=x-xc,yd=y-yc,zd=z-zc,pa=(xd*xa+yd*ya+zd*za)*asi;
        xd-=xa*pa;yd-=ya*pa;zd-=za*pa;
        pa=xd*xd+yd*yd+zd*zd;
        if(pa>1e-5) {
            pa=2*(Math.sqrt(pa)*rc-pa);
            return c.nplane(xd,yd,zd,pa,w_id);
        }
        return true;
    }

    public boolean cut_cell(VoronoiCell  c,double x,double y,double z) {return cut_cell_base(c,x,y,z);}


    public boolean cut_cell(VoronoiCellNeighbor c,double x,double y,double z) {return cut_cell_base(c,x,y,z);}

    private final int w_id;
    private final double xc,yc,zc,xa,ya,za,asi,rc;

}
