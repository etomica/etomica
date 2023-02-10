// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro;


import etomica.util.collections.IntArrayList;

import java.io.OutputStream;
import java.util.Arrays;

import static etomica.util.voro.Common.voro_print_vector;

/** \brief Extension of the voronoicell_base class to represent a Voronoi cell
 * with neighbor information.
 *
 * This class is an extension of the voronoicell_base class, in cases when the
 * IDs of neighboring particles associated with each face of the Voronoi cell.
 * It contains additional data structures mne and ne for storing this
 * information. */
public class VoronoiCellNeighbor extends VoronoiCellBase {
    /** This two dimensional array holds the neighbor information
     * associated with each vertex. mne[p] is a one dimensional
     * array which holds all of the neighbor information for
     * vertices of order p. */
    public int[][] mne;
    /** This is a two dimensional array that holds the neighbor
     * information associated with each vertex. ne[i] points to a
     * one-dimensional array in mne[nu[i]]. ne[i][j] holds the
     * neighbor information associated with the jth edge of vertex
     * i. It is set to the ID number of the plane that made the
     * face that is clockwise from the jth edge. */
//    public int[][] ne;
    public Int2Darray ne_;
    public VoronoiCellNeighbor() {
        super(Config.default_length*Config.default_length);
        memory_setup();
    }
    public VoronoiCellNeighbor(double max_len_sq_) {
        super(max_len_sq_);
        memory_setup();
    }
    public VoronoiCellNeighbor(ContainerBaseBase con) {
        super(con.max_len_sq);
        memory_setup();
    }

    /** Copies the information from another voronoicell class into this
     * class, extending memory allocation if necessary.
     * \param[in] c the class to copy. */
    public void equalOp(VoronoiCell c) {
        check_memory_for_copy(c);
        copy(c);
        int i,j;
        for(i=0;i<c.current_vertex_order;i++) {
            for(j=0;j<c.mec[i]*i;j++) mne[i][j]=0;
            for(j=0;j<c.mec[i];j++) ne_.setIndex(c.mep[i][(2*i+1)*j+2*i], mne[i], j*i);
        }
    }
    public void equalOp(VoronoiCellNeighbor c) {
        check_memory_for_copy(c);
        copy(c);
        int i,j;
        for(i=0;i<c.current_vertex_order;i++) {
            for(j=0;j<c.mec[i]*i;j++) mne[i][j]=c.mne[i][j];
            for(j=0;j<c.mec[i];j++) ne_.setIndex(c.mep[i][(2*i+1)*j+2*i], mne[i], j*i);
        }
    }

    /** This routine calculates the modulus squared of the vector
     * before passing it to the main nplane() routine with full
     * arguments.
     * \param[in] (x,y,z) the vector to cut the cell by.
     * \param[in] p_id the plane ID (for neighbor tracking only).
     * \return False if the plane cut deleted the cell entirely,
     *         true otherwise. */
    public boolean nplane(double x,double y,double z,int p_id) {
        double rsq=x*x+y*y+z*z;
        return nplane(x,y,z,rsq,p_id);
    }
    /** This version of the plane routine just makes up the plane
     * ID to be zero. It will only be referenced if neighbor
     * tracking is enabled.
     * \param[in] (x,y,z) the vector to cut the cell by.
     * \param[in] rsq the modulus squared of the vector.
     * \return False if the plane cut deleted the cell entirely,
     *         true otherwise. */
    public boolean plane(double x,double y,double z,double rsq) {
        return nplane(x,y,z,rsq,0);
    }
    /** Cuts a Voronoi cell using the influence of a particle at
     * (x,y,z), first calculating the modulus squared of this
     * vector before passing it to the main nplane() routine. Zero
     * is supplied as the plane ID, which will be ignored unless
     * neighbor tracking is enabled.
     * \param[in] (x,y,z) the vector to cut the cell by.
     * \return False if the plane cut deleted the cell entirely,
     *         true otherwise. */
    public boolean plane(double x,double y,double z) {
        double rsq=x*x+y*y+z*z;
        return nplane(x,y,z,rsq,0);
    }

    /** This initializes the class to be a rectangular box. It calls the base class
     * initialization routine to set up the edge and vertex information, and then
     * sets up the neighbor information, with initial faces being assigned ID
     * numbers from -1 to -6.
     * \param[in] (xmin,xmax) the minimum and maximum x coordinates.
     * \param[in] (ymin,ymax) the minimum and maximum y coordinates.
     * \param[in] (zmin,zmax) the minimum and maximum z coordinates. */
    public void init(double xmin,double xmax,double ymin,double ymax,double zmin,double zmax) {
        init_base(xmin,xmax,ymin,ymax,zmin,zmax);
        int[] q=mne[3];
    	q[0]=-5;q[1]=-3;q[2]=-1;
        q[3]=-5;q[4]=-2;q[5]=-3;
        q[6]=-5;q[7]=-1;q[8]=-4;
        q[9]=-5;q[10]=-4;q[11]=-2;
        q[12]=-6;q[13]=-1;q[14]=-3;
        q[15]=-6;q[16]=-3;q[17]=-2;
        q[18]=-6;q[19]=-4;q[20]=-1;
        q[21]=-6;q[22]=-2;q[23]=-4;
        ne_.setIndex(0, mne[3], 0);
        ne_.setIndex(1, mne[3], 3);
        ne_.setIndex(2, mne[3], 6);
        ne_.setIndex(3, mne[3], 9);
        ne_.setIndex(4, mne[3], 12);
        ne_.setIndex(5, mne[3], 15);
        ne_.setIndex(6, mne[3], 18);
        ne_.setIndex(7, mne[3], 21);
    }

    /** This initializes the class to be an octahedron. It calls the base class
     * initialization routine to set up the edge and vertex information, and then
     * sets up the neighbor information, with the initial faces being assigned ID
     * numbers from -1 to -8.
     * \param[in] l The distance from the octahedron center to a vertex. Six
     *              vertices are initialized at (-l,0,0), (l,0,0), (0,-l,0),
     *              (0,l,0), (0,0,-l), and (0,0,l). */
    public void init_octahedron(double l) {
        init_octahedron_base(l);
        int[] q=mne[4];
	    q[0]=-5;q[1]=-6;q[2]=-7;q[3]=-8;
        q[4]=-1;q[5]=-2;q[6]=-3;q[7]=-4;
        q[8]=-6;q[9]=-5;q[10]=-2;q[11]=-1;
        q[12]=-8;q[13]=-7;q[14]=-4;q[15]=-3;
        q[16]=-5;q[17]=-8;q[18]=-3;q[19]=-2;
        q[20]=-7;q[21]=-6;q[22]=-1;q[23]=-4;
        ne_.setIndex(0, mne[4], 0);
        ne_.setIndex(1, mne[4], 4);
        ne_.setIndex(2, mne[4], 8);
        ne_.setIndex(3, mne[4], 12);
        ne_.setIndex(4, mne[4], 16);
        ne_.setIndex(5, mne[4], 20);
    }

    /** This initializes the class to be a tetrahedron. It calls the base class
     * initialization routine to set up the edge and vertex information, and then
     * sets up the neighbor information, with the initial faces being assigned ID
     * numbers from -1 to -4.
     * \param (x0,y0,z0) a position vector for the first vertex.
     * \param (x1,y1,z1) a position vector for the second vertex.
     * \param (x2,y2,z2) a position vector for the third vertex.
     * \param (x3,y3,z3) a position vector for the fourth vertex. */
    public void init_tetrahedron(double x0,double y0,double z0,double x1,double y1,double z1,double x2,double y2,double z2,double x3,double y3,double z3) {
        init_tetrahedron_base(x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3);
        int[] q=mne[3];
	    q[0]=-4;q[1]=-3;q[2]=-2;
        q[3]=-3;q[4]=-4;q[5]=-1;
        q[6]=-4;q[7]=-2;q[8]=-1;
        q[9]=-2;q[10]=-3;q[11]=-1;
        ne_.setIndex(0, mne[3], 0);
        ne_.setIndex(1, mne[3], 3);
        ne_.setIndex(2, mne[3], 6);
        ne_.setIndex(3, mne[3], 9);
    }

    /** This routine checks to make sure the neighbor information of each face is
     * consistent. */
    public void check_facets() {
        for(int i=1;i<p;i++) for(int j=0;j<nu[i];j++) {
            int k=ed_.get(i,j);
            if(k>=0) {
                ed_.set(i,j,-1-k);
                int q=ne_.get(i,j);
                int l=cycle_up(ed_.get(i,nu[i]+j), k);
                do {
                    int m=ed_.get(k,l);
                    ed_.set(k,l,-1-m);
                    if(ne_.get(k,l)!=q) System.err.printf("Facet error at (%d,%d)=%d, started from (%d,%d)=%d\n",k,l,ne_.get(k,l),i,j,q);
                    l=cycle_up(ed_.get(k,nu[k]+l),m);
                    k=m;
                } while (k!=i);
            }
        }
        reset_edges();
    }

    /** Computes a vector list of neighbors. */
    public void neighbors(IntArrayList v) {
        v.clear();
        int i,j,k,l,m;
        for(i=1;i<p;i++) for(j=0;j<nu[i];j++) {
            k=ed_.get(i,j);
            if(k>=0) {
                v.add(ne_.get(i,j));
                ed_.set(i,j,-1-k);
                l=cycle_up(ed_.get(i,nu[i]+j),k);
                do {
                    m=ed_.get(k,l);
                    ed_.set(k,l,-1-m);
                    l=cycle_up(ed_.get(k,nu[k]+l),m);
                    k=m;
                } while (k!=i);
            }
        }
        reset_edges();
    }

    /** This prints out the neighbor information for vertex i. */
    public void print_edges_neighbors(int i) {
        if(nu[i]>0) {
            System.out.print("     (");
            for (int j=0; j<nu[i]-1; j++) {
                System.out.printf("%d,",ne_.get(i,j));
            }
            System.out.printf("%d)",ne_.get(i,nu[i]-1));
        } else System.out.print("     ()");
    }
    public void output_neighbors() {
        output_neighbors(System.out);
    }
    public void output_neighbors(OutputStream fp) {
        IntArrayList v = new IntArrayList();
        neighbors(v);
        voro_print_vector(v,fp);
    }

    private Int1Darray paux1_, paux2_;

    /** The class constructor allocates memory for storing neighbor information. */
    private void memory_setup() {
        int i;
        mne=new int[current_vertex_order][];
        ne_ = new Int2Darray();
        for(i=0;i<3;i++) mne[i]=new int[Config.init_n_vertices*i];
        mne[3]=new int[Config.init_3_vertices*3];
        for(i=4;i<current_vertex_order;i++) mne[i]=new int[Config.init_n_vertices*i];
    }
    protected void n_allocate(int i,int m) {mne[i]=new int[m*i];}
    protected void n_add_memory_vertices(int i) {

    }
    protected void n_add_memory_vorder(int i) {
        mne = Arrays.copyOf(mne, i);
    }
    protected void n_set_pointer(int p,int n) {
        ne_.setIndex(p, mne[n], n*mec[n]);
    }
    protected void n_copy(int a,int b,int c,int d) {ne_.set(a,b,ne_.get(c,d));}
    protected void n_set(int a,int b,int c) {ne_.set(a,b,c);}
    protected void n_set_aux1(int k) {
        paux1_ = new Int1Darray(mne[k], k*mec[k]);
    }
    protected void n_copy_aux1(int a,int b) {
        paux1_.set(b, ne_.get(a,b));
    }
    protected void n_copy_aux1_shift(int a,int b) {
        paux1_.set(b, ne_.get(a, b+1));
    }
    protected void n_set_aux2_copy(int a,int b) {
        paux2_ = new Int1Darray(mne[b], b*mec[b]);
        for(int i=0;i<b;i++) ne_.set(a,i,paux2_.get(i));
    }
    protected void n_copy_pointer(int a,int b) {ne_.copyOffset(b,a);}
    protected void n_set_to_aux1(int j) {ne_.setIndex(j, paux1_.storage, paux1_.getOffset());}
    protected void n_set_to_aux2(int j) {ne_.setIndex(j, paux2_.storage, paux2_.offset);}
    protected void n_allocate_aux1(int i) {paux1_=new Int1Darray(new int[i*mem[i]],0);}
    protected void n_switch_to_aux1(int i) {
        if (paux1_.offset!=0) throw new RuntimeException("oops");
        mne[i]=paux1_.storage;
    }
    protected void n_copy_to_aux1(int i,int m) {paux1_.set(i, mne[i][m]);}
    protected void n_set_to_aux1_offset(int k,int m) {ne_.setIndex(k, paux1_.storage, paux1_.getOffset()+m);}

}
