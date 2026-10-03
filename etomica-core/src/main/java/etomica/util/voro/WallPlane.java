// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro;

/** \brief A class representing a plane wall object.
 *
 * This class represents a single plane wall object. */
public class WallPlane implements Wall {

    public WallPlane(double xc_, double yc_, double zc_, double ac_) {
        this(xc_, yc_, zc_, ac_, -99);
    }

    /** Constructs a plane wall object.
     * \param[in] (xc_,yc_,zc_) a normal vector to the plane.
     * \param[in] ac_ a displacement along the normal vector.
     * \param[in] w_id_ an ID number to associate with the wall for
     *		    neighbor tracking. */
    public WallPlane(double xc_,double yc_,double zc_,double ac_,int w_id_) {
        w_id = w_id_;
        xc = xc_;
        yc = yc_;
        zc = zc_;
        ac = ac_;
    }
    public boolean point_inside(double x,double y,double z) {
        return x*xc+y*yc+z*zc<ac;
    }

    public boolean cut_cell_base(VoronoiCellBase c,double x,double y,double z) {
        double dq=2*(ac-x*xc-y*yc-z*zc);
        return c.nplane(xc,yc,zc,dq,w_id);
    }
    public boolean cut_cell(VoronoiCell c,double x,double y,double z) {return cut_cell_base(c,x,y,z);}
    public boolean cut_cell(VoronoiCellNeighbor c,double x,double y,double z) {return cut_cell_base(c,x,y,z);}

    private final int w_id;
    private final double xc,yc,zc,ac;
}
