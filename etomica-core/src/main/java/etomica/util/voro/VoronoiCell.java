// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro;

import etomica.util.collections.IntArrayList;

/** \brief Extension of the voronoicell_base class to represent a Voronoi
 * cell without neighbor information.
 *
 * This class is an extension of the voronoicell_base class, in cases when
 * is not necessary to track the IDs of neighboring particles associated
 * with each face of the Voronoi cell. */
public class VoronoiCell extends VoronoiCellBase {

    public VoronoiCell() {
        super(Config.default_length*Config.default_length);
    }

    public VoronoiCell(double max_len_sq_) {
        super(max_len_sq_);
    }

    public VoronoiCell(ContainerBaseBase con) {
        this(con.max_len_sq);
    }

    public VoronoiCell(VoronoiCell c) {
        this();
        equalOperator(c);
    }

    /** Copies the information from another voronoicell class into
     * this class, extending memory allocation if necessary.
     * \param[in] c the class to copy. */
    public void equalOperator(VoronoiCell c) {
        check_memory_for_copy(c);
        copy(c);
    }
    /** Cuts a Voronoi cell using by the plane corresponding to the
     * perpendicular bisector of a particle.
     * \param[in] (x,y,z) the position of the particle.
     * \param[in] rsq the modulus squared of the vector.
     * \param[in] p_id the plane ID, ignored for this case where no
     *                 neighbor tracking is enabled.
     * \return False if the plane cut deleted the cell entirely,
     *         true otherwise. */
    public boolean nplane(double x,double y,double z,double rsq,int p_id) {
        return super.nplane(x,y,z,rsq,0);
    }
    /** Cuts a Voronoi cell using by the plane corresponding to the
     * perpendicular bisector of a particle.
     * \param[in] (x,y,z) the position of the particle.
     * \param[in] p_id the plane ID, ignored for this case where no
     *                 neighbor tracking is enabled.
     * \return False if the plane cut deleted the cell entirely,
     *         true otherwise. */
    public boolean nplane(double x,double y,double z,int p_id) {
        double rsq=x*x+y*y+z*z;
        return super.nplane(x,y,z,rsq,0);
    }
    /** Cuts a Voronoi cell using by the plane corresponding to the
     * perpendicular bisector of a particle.
     * \param[in] (x,y,z) the position of the particle.
     * \param[in] rsq the modulus squared of the vector.
     * \return False if the plane cut deleted the cell entirely,
     *         true otherwise. */
    public boolean plane(double x,double y,double z,double rsq) {
        return super.nplane(x,y,z,rsq,0);
    }
    /** Cuts a Voronoi cell using by the plane corresponding to the
     * perpendicular bisector of a particle.
     * \param[in] (x,y,z) the position of the particle.
     * \return False if the plane cut deleted the cell entirely,
     *         true otherwise. */
    public boolean plane(double x,double y,double z) {
        double rsq=x*x+y*y+z*z;
        return nplane(x,y,z,rsq,0);
    }
    /** Initializes the Voronoi cell to be rectangular box with the
     * given dimensions.
     * \param[in] (xmin,xmax) the minimum and maximum x coordinates.
     * \param[in] (ymin,ymax) the minimum and maximum y coordinates.
     * \param[in] (zmin,zmax) the minimum and maximum z coordinates. */
    public void init(double xmin,double xmax,double ymin,double ymax,double zmin,double zmax) {
        init_base(xmin,xmax,ymin,ymax,zmin,zmax);
    }
    /** Initializes the cell to be an octahedron with vertices at
     * (l,0,0), (-l,0,0), (0,l,0), (0,-l,0), (0,0,l), and (0,0,-l).
     * \param[in] l a parameter setting the size of the octahedron.
     */
    public void init_octahedron(double l) {
        init_octahedron_base(l);
    }
    /** Initializes the cell to be a tetrahedron.
     * \param[in] (x0,y0,z0) the coordinates of the first vertex.
     * \param[in] (x1,y1,z1) the coordinates of the second vertex.
     * \param[in] (x2,y2,z2) the coordinates of the third vertex.
     * \param[in] (x3,y3,z3) the coordinates of the fourth vertex.
     */
    public void init_tetrahedron(double x0,double y0,double z0,double x1,double y1,double z1,double x2,double y2,double z2,double x3,double y3,double z3) {
        init_tetrahedron_base(x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3);
    }

    public void init_l_shape() {
        for(int i=0;i<current_vertex_order;i++) mec[i]=0;
        up=0;
        mec[3]=p=12;
    	final double j=0;
	    pts[0]=-2;pts[1]=-2;pts[2]=-2;
        pts[4]=2;pts[5]=-2;pts[6]=-2;
        pts[8]=-2;pts[9]=0;pts[10]=-2;
        pts[12]=-j;pts[13]=j;pts[14]=-2;
        pts[16]=0;pts[17]=2;pts[18]=-2;
        pts[20]=2;pts[21]=2;pts[22]=-2;
        pts[24]=-2;pts[25]=-2;pts[26]=2;
        pts[28]=2;pts[29]=-2;pts[30]=2;
        pts[32]=-2;pts[33]=0;pts[34]=2;
        pts[36]=-j;pts[37]=j;pts[38]=2;
        pts[40]=0;pts[41]=2;pts[42]=2;
        pts[44]=2;pts[45]=2;pts[46]=2;
        int[] q=mep[3];
	    q[0]=1;q[1]=6;q[2]=2;q[6]=0;
        q[7]=5;q[8]=7;q[9]=0;q[13]=1;
        q[14]=0;q[15]=8;q[16]=3;q[20]=2;
        q[21]=2;q[22]=9;q[23]=4;q[27]=3;
        q[28]=3;q[29]=10;q[30]=5;q[34]=4;
        q[35]=4;q[36]=11;q[37]=1;q[41]=5;
        q[42]=8;q[43]=0;q[44]=7;q[48]=6;
        q[49]=6;q[50]=1;q[51]=11;q[55]=7;
        q[56]=9;q[57]=2;q[58]=6;q[62]=8;
        q[63]=10;q[64]=3;q[65]=8;q[69]=9;
        q[70]=11;q[71]=4;q[72]=9;q[76]=10;
        q[77]=7;q[78]=5;q[79]=10;q[83]=11;
        ed_.setOffsets(q, new int[]{0,7,14,21,28,35,42,49,56,63,70,77});
        for(int i=0;i<12;i++) nu[i]=3;
        construct_relations();
    }

    protected void n_allocate(int i,int m) {};
    protected void n_add_memory_vertices(int i) {};
    protected void n_add_memory_vorder(int i) {};
    protected void n_set_pointer(int p,int n) {};
    protected void n_copy(int a,int b,int c,int d) {};
    protected void n_set(int a,int b,int c) {};
    protected void n_set_aux1(int k) {};
    protected void n_copy_aux1(int a,int b) {};
    protected void n_copy_aux1_shift(int a,int b) {};
    protected void n_set_aux2_copy(int a,int b) {};
    protected void n_copy_pointer(int a,int b) {};
    protected void n_set_to_aux1(int j) {};
    protected void n_set_to_aux2(int j) {};
    protected void n_allocate_aux1(int i) {};
    protected void n_switch_to_aux1(int i) {};
    protected void n_copy_to_aux1(int i,int m) {};
    protected void n_set_to_aux1_offset(int k,int m) {};
    protected void n_neighbors(IntArrayList v) {v.clear();};

}
