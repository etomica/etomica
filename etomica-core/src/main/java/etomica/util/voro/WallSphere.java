// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro;

/** \brief A class representing a spherical wall object.
 *
 * This class represents a spherical wall object. */
public class WallSphere implements Wall {

    public WallSphere(double xc_, double yc_, double zc_, double rc_) {
        this(xc_, yc_, zc_, rc_, -99);
    }

    /** Constructs a spherical wall object.
     * \param[in] w_id_ an ID number to associate with the wall for
     *		    neighbor tracking.
     * \param[in] (xc_,yc_,zc_) a position vector for the sphere's
     * 			    center.
     * \param[in] rc_ the radius of the sphere. */
    public WallSphere(double xc_,double yc_,double zc_,double rc_,int w_id_) {
        w_id = w_id_;
        xc = xc_;
        yc = yc_;
        zc = zc_;
        rc = rc_;
    }

    /** Tests to see whether a point is inside the sphere wall object.
     * \param[in,out] (x,y,z) the vector to test.
     * \return True if the point is inside, false if the point is outside. */
    public boolean point_inside(double x,double y,double z){
        return (x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc)<rc*rc;
    }

    /** Cuts a cell by the sphere wall object. The spherical wall is approximated by
     * a single plane applied at the point on the sphere which is closest to the center
     * of the cell. This works well for particle arrangements that are packed against
     * the wall, but loses accuracy for sparse particle distributions.
     * \param[in,out] c the Voronoi cell to be cut.
     * \param[in] (x,y,z) the location of the Voronoi cell.
     * \return True if the cell still exists, false if the cell is deleted. */
    public boolean cut_cell_base(VoronoiCellBase c,double x,double y,double z) {
        double xd=x-xc,yd=y-yc,zd=z-zc,dq=xd*xd+yd*yd+zd*zd;
        if (dq>1e-5) {
            dq=2*(Math.sqrt(dq)*rc-dq);
            return c.nplane(xd,yd,zd,dq,w_id);
        }
        return true;
    }
    public boolean cut_cell(VoronoiCell c,double x,double y,double z) {return cut_cell_base(c,x,y,z);}
    public boolean cut_cell(VoronoiCellNeighbor c,double x,double y,double z) {return cut_cell_base(c,x,y,z);}

    private final int w_id;
    private double xc,yc,zc,rc;
}
