// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro;

import etomica.exception.MethodNotImplementedException;
import etomica.util.collections.DoubleArrayList;
import etomica.util.collections.IntArrayList;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Arrays;

import static etomica.util.voro.Common.*;

/** \brief A class representing a single Voronoi cell.
 *
 * This class represents a single Voronoi cell, as a collection of vertices
 * that are connected by edges. The class contains routines for initializing
 * the Voronoi cell to be simple shapes such as a box, tetrahedron, or octahedron.
 * It the contains routines for recomputing the cell based on cutting it
 * by a plane, which forms the key routine for the Voronoi cell computation.
 * It contains numerous routine for computing statistics about the Voronoi cell,
 * and it can output the cell in several formats.
 *
 * This class is not intended for direct use, but forms the base of the
 * voronoicell and voronoicell_neighbor classes, which extend it based on
 * whether neighboring particle ID information needs to be tracked. */
public abstract class VoronoiCellBase {

    /** This holds the current size of the arrays ed and nu, which
     * hold the vertex information. If more vertices are created
     * than can fit in this array, then it is dynamically extended
     * using the add_memory_vertices routine. */
    public int current_vertices;
    /** This holds the current maximum allowed order of a vertex,
     * which sets the size of the mem, mep, and mec arrays. If a
     * vertex is created with more vertices than this, the arrays
     * are dynamically extended using the add_memory_vorder routine.
     */
    public int current_vertex_order;
    /** This sets the size of the main delete stack. */
    public int current_delete_size;
    /** This sets the size of the auxiliary delete stack. */
    public int current_delete2_size;
    /** This sets the size of the extra search stack. */
    int current_xsearch_size;
    /** This sets the total number of vertices in the current cell.
     */
    public int p;
    /** This is the index of particular point in the cell, which is
     * used to start the tracing routines for plane intersection
     * and cutting. These routines will work starting from any
     * point, but it's often most efficient to start from the last
     * point considered, since in many cases, the cell construction
     * algorithm may consider many planes with similar vectors
     * concurrently. */
    public int up;
    /** This is a two dimensional array that holds information
     * about the edge connections of the vertices that make up the
     * cell. The two dimensional array is not allocated in the
     * usual method. To account for the fact the different vertices
     * have different orders, and thus require different amounts of
     * storage, the elements of ed[i] point to one-dimensional
     * arrays in the mep[] array of different sizes.
     *
     * More specifically, if vertex i has order m, then ed[i]
     * points to a one-dimensional array in mep[m] that has 2*m+1
     * entries. The first m elements hold the neighboring edges, so
     * that the jth edge of vertex i is held in ed[i][j]. The next
     * m elements hold a table of relations which is redundant but
     * helps speed up the computation. It satisfies the relation
     * ed[ed[i][j]][ed[i][m+j]]=i. The final entry holds a back
     * pointer, so that ed[i+2*m]=i. The back pointers are used
     * when rearranging the memory. */
//    public int[][] ed;
    public Int2Darray ed_;
    /** This array holds the order of the vertices in the Voronoi
     * cell. This array is dynamically allocated, with its current
     * size held by current_vertices. */
    public int[] nu;
    public int[] mask;
    /** This in an array with size 3*current_vertices for holding
     * the positions of the vertices. */
    public double[] pts;
    public double tol;
    public double tol_cu;
    public double big_tol;

    /** Constructs a Voronoi cell and sets up the initial memory. */
    public VoronoiCellBase(double max_len_sq) {
        current_vertices = Config.init_vertices;
        current_vertex_order = Config.init_vertex_order;
        current_delete_size = Config.init_delete_size;
        current_delete2_size = Config.init_delete2_size;
        current_xsearch_size = Config.init_xsearch_size;
        nu = new int[current_vertices];
        mask = new int[current_vertices];
        pts = new double[4*current_vertices];
        tol = Config.tolerance*max_len_sq;
        tol_cu = tol*Math.sqrt(tol);
        big_tol = Config.big_tolerance_fac*tol;
        mem = new int[current_vertex_order];
        mec = new int[current_vertex_order];
        mep = new int[current_vertex_order][];
        ds = new int[current_delete_size];
        ds2 = new int[current_delete2_size];
        xse = new int[current_xsearch_size];
        maskc = 0;
        for(int i=0;i<3;i++) {
            mem[i]=Config.init_n_vertices;mec[i]=0;
            mep[i]=new int[Config.init_n_vertices*((i<<1)+1)];
        }
        mem[3]=Config.init_3_vertices;
        mec[3]=0;
        mep[3]=new int[Config.init_3_vertices*7];
        for(int i=4;i<current_vertex_order;i++) {
            mem[i]=Config.init_n_vertices;
            mec[i]=0;
            mep[i]=new int[Config.init_n_vertices*((i<<1)+1)];
        }

        ed_ = new Int2Darray();

    }
    public abstract void init(double xmin,double xmax,double ymin,double ymax,double zmin,double zmax);

    /** Initializes a Voronoi cell as a rectangular box with the given dimensions.
     * \param[in] (xmin,xmax) the minimum and maximum x coordinates.
     * \param[in] (ymin,ymax) the minimum and maximum y coordinates.
     * \param[in] (zmin,zmax) the minimum and maximum z coordinates. */
    public void init_base(double xmin,double xmax,double ymin,double ymax,double zmin,double zmax) {
        for(int i=0;i<current_vertex_order;i++) mec[i]=0;
        up=0;
        mec[3]=p=8;xmin*=2;xmax*=2;ymin*=2;ymax*=2;zmin*=2;zmax*=2;
	    pts[0]=xmin;pts[1]=ymin;pts[2]=zmin;
        pts[4]=xmax;pts[5]=ymin;pts[6]=zmin;
        pts[8]=xmin;pts[9]=ymax;pts[10]=zmin;
        pts[12]=xmax;pts[13]=ymax;pts[14]=zmin;
        pts[16]=xmin;pts[17]=ymin;pts[18]=zmax;
        pts[20]=xmax;pts[21]=ymin;pts[22]=zmax;
        pts[24]=xmin;pts[25]=ymax;pts[26]=zmax;
        pts[28]=xmax;pts[29]=ymax;pts[30]=zmax;
        int[] q=mep[3];
        q[0]=1;q[1]=4;q[2]=2;q[3]=2;q[4]=1;q[5]=0;q[6]=0;
        q[7]=3;q[8]=5;q[9]=0;q[10]=2;q[11]=1;q[12]=0;q[13]=1;
        q[14]=0;q[15]=6;q[16]=3;q[17]=2;q[18]=1;q[19]=0;q[20]=2;
        q[21]=2;q[22]=7;q[23]=1;q[24]=2;q[25]=1;q[26]=0;q[27]=3;
        q[28]=6;q[29]=0;q[30]=5;q[31]=2;q[32]=1;q[33]=0;q[34]=4;
        q[35]=4;q[36]=1;q[37]=7;q[38]=2;q[39]=1;q[40]=0;q[41]=5;
        q[42]=7;q[43]=2;q[44]=4;q[45]=2;q[46]=1;q[47]=0;q[48]=6;
        q[49]=5;q[50]=3;q[51]=6;q[52]=2;q[53]=1;q[54]=0;q[55]=7;
        ed_.setOffsets(mep[3], new int[]{0,7,14,21,28,35,42,49});
	    nu[0]=nu[1]=nu[2]=nu[3]=nu[4]=nu[5]=nu[6]=nu[7]=3;
    }

    /** Initializes an L-shaped Voronoi cell of a fixed size for testing the
     * convexity robustness. */
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
         ed_.setOffsets(mep[3], new int[]{0,7,14,21,28,35,42,49,56,63,70,77});
         for(int i=0;i<12;i++) nu[i]=3;
         construct_relations();
     }

    /** Initializes a Voronoi cell as a regular octahedron.
     * \param[in] l The distance from the octahedron center to a vertex. Six
     *              vertices are initialized at (-l,0,0), (l,0,0), (0,-l,0),
     *              (0,l,0), (0,0,-l), and (0,0,l). */
    public void init_octahedron_base(double l) {
        for(int i=0;i<current_vertex_order;i++) mec[i]=0;
        up=0;
        mec[4]=p=6;l*=2;
        pts[0]=-l;pts[1]=0;pts[2]=0;
        pts[4]=l;pts[5]=0;pts[6]=0;
        pts[8]=0;pts[9]=-l;pts[10]=0;
        pts[12]=0;pts[13]=l;pts[14]=0;
        pts[16]=0;pts[17]=0;pts[18]=-l;
        pts[20]=0;pts[21]=0;pts[22]=l;
        int[] q=mep[4];
	    q[0]=2;q[1]=5;q[2]=3;q[3]=4;q[4]=0;q[5]=0;q[6]=0;q[7]=0;q[8]=0;
        q[9]=2;q[10]=4;q[11]=3;q[12]=5;q[13]=2;q[14]=2;q[15]=2;q[16]=2;q[17]=1;
        q[18]=0;q[19]=4;q[20]=1;q[21]=5;q[22]=0;q[23]=3;q[24]=0;q[25]=1;q[26]=2;
        q[27]=0;q[28]=5;q[29]=1;q[30]=4;q[31]=2;q[32]=3;q[33]=2;q[34]=1;q[35]=3;
        q[36]=0;q[37]=3;q[38]=1;q[39]=2;q[40]=3;q[41]=3;q[42]=1;q[43]=1;q[44]=4;
        q[45]=0;q[46]=2;q[47]=1;q[48]=3;q[49]=1;q[50]=3;q[51]=3;q[52]=1;q[53]=5;
        ed_.setOffsets(mep[4], new int[]{0,9,18,27,36,45});
//	    ed[0]=q;ed[1]=q+9;ed[2]=q+18;ed[3]=q+27;ed[4]=q+36;ed[5]=q+45;
	    nu[0]=nu[1]=nu[2]=nu[3]=nu[4]=nu[5]=4;
    }

    /** Initializes a Voronoi cell as a tetrahedron. It assumes that the normal to
     * the face for the first three vertices points inside.
     * \param (x0,y0,z0) a position vector for the first vertex.
     * \param (x1,y1,z1) a position vector for the second vertex.
     * \param (x2,y2,z2) a position vector for the third vertex.
     * \param (x3,y3,z3) a position vector for the fourth vertex. */
    public void init_tetrahedron_base(double x0,double y0,double z0,double x1,double y1,double z1,double x2,double y2,double z2,double x3,double y3,double z3) {
        for(int i=0;i<current_vertex_order;i++) mec[i]=0;
        up=0;
        mec[3]=p=4;
    	pts[0]=x0*2;pts[1]=y0*2;pts[2]=z0*2;
        pts[3]=x1*2;pts[4]=y1*2;pts[5]=z1*2;
        pts[6]=x2*2;pts[7]=y2*2;pts[8]=z2*2;
        pts[9]=x3*2;pts[10]=y3*2;pts[11]=z3*2;
        int[] q=mep[3];
	    q[0]=1;q[1]=3;q[2]=2;q[3]=0;q[4]=0;q[5]=0;q[6]=0;
        q[7]=0;q[8]=2;q[9]=3;q[10]=0;q[11]=2;q[12]=1;q[13]=1;
        q[14]=0;q[15]=3;q[16]=1;q[17]=2;q[18]=2;q[19]=1;q[20]=2;
        q[21]=0;q[22]=1;q[23]=2;q[24]=1;q[25]=2;q[26]=1;q[27]=3;
        ed_.setOffsets(mep[3], new int[]{0,7,14,21});
	    nu[0]=nu[1]=nu[2]=nu[3]=3;
    }

    /** Translates the vertices of the Voronoi cell by a given vector.
     * \param[in] (x,y,z) the coordinates of the vector. */
    public void translate(double x,double y,double z) {
        x*=2;y*=2;z*=2;
        for (int i=0; i<p; i+=4) {
            pts[i+0] = x;
            pts[i+1] = y;
            pts[i+2] = z;
        }
    }

    /** Outputs the edges of the Voronoi cell in POV-Ray format to an open file
     * stream, displacing the cell by given vector.
     * \param[in] (x,y,z) a displacement vector to be added to the cell's position.
     * \param[in] fp a file handle to write to. */
    public void draw_pov(double x, double y, double z) {
        draw_pov(x, y, z, System.out);
    }
    public void draw_pov(double x, double y, double z, OutputStream fp) {
        try {
            for (int i = 0; i < p; i++) {
                String posbuf1 = String.format("%g,%g,%g", x + pts[4 * i + 0] * 0.5, y + pts[4 * i + 1] * 0.5, z + pts[4 * i + 2] * 0.5);
                fp.write(String.format("sphere{<%s>,r}\n", posbuf1).getBytes());
                for (int j = 0; j < nu[i]; j++) {
                    int k = ed_.get(i,j);
                    if (k < i) {
                        String posbuf2 = String.format("%g,%g,%g", x + pts[4 * k + 0] * 0.5, y + 0.5 * pts[4 * k + 1], z + 0.5 * pts[4 * k + 2]);
                        if (!posbuf1.equals(posbuf2))
                            fp.write(("cylinder{<" + posbuf1 + ">,<" + posbuf2 + ">,r}\n").getBytes());
                    }
                }
            }
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }
    /** Outputs the cell in POV-Ray format, using cylinders for edges
     * and spheres for vertices, to a given file.
     * \param[in] (x,y,z) a displacement to add to the cell's
     *                    position.
     * \param[in] filename the name of the file to write to. */
    public void draw_pov(double x,double y,double z, String filename) {
        try {
            FileOutputStream fw = new FileOutputStream(filename);
            draw_pov(x, y, z, fw);
            fw.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    };
    public void draw_pov_mesh(double x, double y, double z) {
        draw_pov_mesh(x, y, z, System.out);
    }
    public void draw_pov_mesh(double x,double y,double z, OutputStream fp) {
        try {
            fp.write(("mesh2 {\nvertex_vectors {\n" + p + "\n").getBytes());
            for (int i = 0; i < p; i++)
                fp.write((",<" + (x + pts[4 * i + 0]) + "," + (y + pts[4 * i + 1]) + "," + (z + pts[4 * i + 2]) + ">\n").getBytes());
            fp.write(("}\nface_indices {\n"+((p - 2) << 1)+"\n").getBytes());
            for (int i = 1; i < p; i++)
                for (int j = 0; j < nu[i]; j++) {
                    int k = ed_.get(i,j);
                    if (k >= 0) {
                        ed_.set(i,j,-1-k);
                        int l = cycle_up(ed_.get(i,nu[i]), k);
                        int m = ed_.get(k,l);
                        ed_.set(k,l,-1-m);
                        while (m != i) {
                            int n = cycle_up(ed_.get(k,nu[k]+l), m);
                            fp.write(String.format(",<%d,%d,%d>\n", i, k, m).getBytes());
                            k = m;
                            l = n;
                            m = ed_.get(k,l);
                            ed_.set(k,l,-1-m);
                        }
                    }
                }
            fp.write("}\ninside_vector <0,0,1>\n}\n".getBytes());
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
        reset_edges();
    }
    /** Outputs the cell in POV-Ray format as a mesh2 object to a
     * given file.
     * \param[in] (x,y,z) a displacement to add to the cell's
     *                    position.
     * \param[in] filename the name of the file to write to. */
    public void draw_pov_mesh(double x,double y,double z, String filename) {
        try {
            FileOutputStream fw = new FileOutputStream(filename);
            draw_pov_mesh(x, y, z, fw);
            fw.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    /** Outputs the edges of the Voronoi cell in gnuplot format to an output stream.
     * \param[in] (x,y,z) a displacement vector to be added to the cell's position.
     * \param[in] fp a file handle to write to. */
    public void draw_gnuplot(double x,double y,double z) {
        draw_gnuplot(x, y, z, System.out);
    }
    public void draw_gnuplot(double x,double y,double z,OutputStream fp) {
        try {
            for (int i = 1; i < p; i++) {
                for (int j = 0; j < nu[i]; j++) {
                    int k = ed_.get(i,j);
                    if (k >= 0) {
                        fp.write(((x + 0.5 * pts[4*i]) + " " + (y + 0.5 * pts[4 * i + 1]) + " " + (z + 0.5 * pts[4 * i + 2]) + "\n").getBytes());
                        int l = i;
                        int m = j;
                        do {
                            ed_.set(k,ed_.get(l,nu[l]+m), -1-l);
                            ed_.set(l,m,-1-k);
                            l = k;
                            fp.write(((x + 0.5 * pts[4 * k]) + " " + (y + 0.5 * pts[4 * k + 1]) + " " + (z + 0.5 * pts[4 * k + 2]) + "\n").getBytes());
                            m = search_edge(l);
                            if (m>=0) k = ed_.get(l,m);
                        } while (m>=0);
                        fp.write("\n\n".getBytes());
                    }
                }
            }
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
        reset_edges();
    }
    /** Outputs the cell in Gnuplot format a given file.
     * \param[in] (x,y,z) a displacement to add to the cell's
     *                    position.
     * \param[in] filename the name of the file to write to. */
    public void draw_gnuplot(double x,double y,double z,String filename) {
        try {
            FileOutputStream fw = new FileOutputStream(filename);
            draw_gnuplot(x, y, z, fw);
            fw.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    /** Calculates the volume of the Voronoi cell, by decomposing the cell into
     * tetrahedra extending outward from the zeroth vertex, whose volumes are
     * evaluated using a scalar triple product.
     * \return A floating point number holding the calculated volume. */
    public double volume() {
	    final double fe=1/48.0;
        double vol=0;
        int i,j,k,l,m,n;
        double ux,uy,uz,vx,vy,vz,wx,wy,wz;
        for(i=1;i<p;i++) {
            ux=pts[0]-pts[i<<2];
            uy=pts[1]-pts[(i<<2)+1];
            uz=pts[2]-pts[(i<<2)+2];
            for(j=0;j<nu[i];j++) {
                k = ed_.get(i,j);
                if(k>=0) {
                    ed_.set(i,j,-1-k);
                    l = cycle_up(ed_.get(i,nu[i]+j), k);
                    vx=pts[k<<2]-pts[0];
                    vy=pts[(k<<2)+1]-pts[1];
                    vz=pts[(k<<2)+2]-pts[2];
                    m = ed_.get(k,l);
                    ed_.set(k,l,-1-m);

                    while(m!=i) {
                        n = cycle_up(ed_.get(k, nu[k]+l), m);
                        wx=pts[m<<2]-pts[0];
                        wy=pts[(m<<2)+1]-pts[1];
                        wz=pts[(m<<2)+2]-pts[2];
                        vol+=ux*vy*wz+uy*vz*wx+uz*vx*wy-uz*vy*wx-uy*vx*wz-ux*vz*wy;
                        k=m;l=n;vx=wx;vy=wy;vz=wz;
                        m = ed_.get(k,l);
                        ed_.set(k,l,-1-m);
                    }
                }
            }
        }
        reset_edges();
        return vol*fe;
    }

    /** Computes the maximum radius squared of a vertex from the center of the
     * cell. It can be used to determine when enough particles have been testing an
     * all planes that could cut the cell have been considered.
     * \return The maximum radius squared of a vertex.*/
    public double max_radius_squared() {
        double r = pts[0]*pts[0] + pts[1]*pts[1] + pts[2]*pts[2];
        for (int i=4; i<4*p; i+=4) {
            double s = pts[i]*pts[i];
            s += pts[i+1]*pts[i+1];
            s += pts[i+2]*pts[i+2];
            r = Math.max(r, s);
        }
        return r;
    }

    /** Calculates the total edge distance of the Voronoi cell.
     * \return A floating point number holding the calculated distance. */
    public double total_edge_distance() {
        int i,j,k;
        double dis=0,dx,dy,dz;
        for(i=0;i<p-1;i++) for(j=0;j<nu[i];j++) {
            k = ed_.get(i,j);
            if(k>i) {
                dx=pts[k<<2]-pts[i<<2];
                dy=pts[(k<<2)+1]-pts[(i<<2)+1];
                dz=pts[(k<<2)+2]-pts[(i<<2)+2];
                dis+=Math.sqrt(dx*dx+dy*dy+dz*dz);
            }
        }
        return 0.5*dis;
    }

    /** Calculates the total surface area of the Voronoi cell.
     * \return The computed area. */
    public double surface_area() {
        double area=0;
        int i,j,k,l,m,n;
        double ux,uy,uz,vx,vy,vz,wx,wy,wz;
        for(i=1;i<p;i++) for(j=0;j<nu[i];j++) {
            k = ed_.get(i,j);
            if(k>=0) {
                ed_.set(i,j,-1-k);
                l = cycle_up(ed_.get(i,nu[i]+j), k);
                m = ed_.get(k,l);
                ed_.set(k,l,-1-m);
                while(m!=i) {
                    n = cycle_up(ed_.get(k,nu[k]+l), m);
                    ux=pts[4*k]-pts[4*i];
                    uy=pts[4*k+1]-pts[4*i+1];
                    uz=pts[4*k+2]-pts[4*i+2];
                    vx=pts[4*m]-pts[4*i];
                    vy=pts[4*m+1]-pts[4*i+1];
                    vz=pts[4*m+2]-pts[4*i+2];
                    wx=uy*vz-uz*vy;
                    wy=uz*vx-ux*vz;
                    wz=ux*vy-uy*vx;
                    area+=Math.sqrt(wx*wx+wy*wy+wz*wz);
                    k=m;l=n;
//                    m=ed_.get(k,l);ed_.set(k,l,-1-m);
                    m = ed_.get(k,l);
                    ed_.set(k,l,-1-m);
                }
            }
        }
        reset_edges();
        return 0.125*area;
    }

    public double[] centroid() {
        double tvol,vol=0;
        double[] cxyz = new double[3];
        int i,j,k,l,m,n;
        double ux,uy,uz,vx,vy,vz,wx,wy,wz;
        for(i=1;i<p;i++) {
            ux=pts[0]-pts[4*i];
            uy=pts[1]-pts[4*i+1];
            uz=pts[2]-pts[4*i+2];
            for(j=0;j<nu[i];j++) {
                k = ed_.get(i,j);
                if(k>=0) {
                    ed_.set(i,j,-1-k);
                    l = cycle_up(ed_.get(i,nu[i]+j),k);
                    vx=pts[4*k]-pts[0];
                    vy=pts[4*k+1]-pts[1];
                    vz=pts[4*k+2]-pts[2];
                    m = ed_.get(k,l);
                    ed_.set(k,l,-1-m);
                    while(m!=i) {
                        n = cycle_up(ed_.get(k,nu[k]+l), m);
                        wx=pts[4*m]-pts[0];
                        wy=pts[4*m+1]-pts[1];
                        wz=pts[4*m+2]-pts[2];
                        tvol=ux*vy*wz+uy*vz*wx+uz*vx*wy-uz*vy*wx-uy*vx*wz-ux*vz*wy;
                        vol+=tvol;
                        cxyz[0]+=(wx+vx-ux)*tvol;
                        cxyz[1]+=(wy+vy-uy)*tvol;
                        cxyz[2]+=(wz+vz-uz)*tvol;
                        k=m;l=n;vx=wx;vy=wy;vz=wz;
                        m = ed_.get(k,l);
                        ed_.set(k,l,-1-m);
                    }
                }
            }
        }
        reset_edges();
        if(vol>tol_cu) {
            vol=0.125/vol;
            cxyz[0]=cxyz[0]*vol+0.5*pts[0];
            cxyz[1]=cxyz[1]*vol+0.5*pts[1];
            cxyz[2]=cxyz[2]*vol+0.5*pts[2];
        } else cxyz[0]=cxyz[1]=cxyz[2]=0;
        return cxyz;
    }

    /** Returns the number of faces of a computed Voronoi cell.
     * \return The number of faces. */
    public int number_of_faces() {
        int i,j,k,l,m,s=0;
        for(i=1;i<p;i++) for(j=0;j<nu[i];j++) {
            k = ed_.get(i,j);
            if(k>=0) {
                s++;
                ed_.set(i,j,-1-k);
                l = cycle_up(ed_.get(i,nu[i]+j), k);
                do {
                    m = ed_.get(k,l);
                    ed_.set(k,l,-1-m);
                    l = cycle_up(ed_.get(k,nu[k]+l), m);
                    k=m;
                } while (k!=i);

            }
        }
        reset_edges();
        return s;
    }

    /** Counts the number of edges of the Voronoi cell.
     * \return the number of edges. */
    public int number_of_edges() {
        int edges=0;
        for (int i=0; i<p; i++) edges += nu[i];
        return edges>>1;
    }

    /** Returns a vector of the vertex orders.
     * \param[out] v the vector to store the results in. */
    public void vertex_orders(IntArrayList v) {
        v.clear();
        for(int i=0;i<p;i++) v.add(nu[i]);
    }

    /** Outputs the vertex orders.
     * \param[out] fp the file handle to write to. */
    public void output_vertex_orders() {
        output_vertex_orders(System.out);
    }
    public void output_vertex_orders(OutputStream fp) {
        if(p>0) {
            try {
                fp.write(("" + nu[0]).getBytes());
                for (int i = 1; i < p; i++) fp.write(("" + nu[0]).getBytes());
            }
            catch (IOException ex) {
                throw new RuntimeException(ex);
            }
        }
    }

    /** Returns a vector of the vertex vectors using the local coordinate system.
     * \param[out] v the vector to store the results in. */
    public void vertices(DoubleArrayList v) {
        v.clear();
        for(int i=0;i<p;i++) {
            v.add(pts[4*i+0]*0.5);
            v.add(pts[4*i+1]*0.5);
            v.add(pts[4*i+2]*0.5);
        }
    }

    /** Outputs the vertex vectors using the local coordinate system.
     * \param[out] fp the file handle to write to. */
    public void output_vertices() {
        output_vertices(System.out);
    }
    public void output_vertices(OutputStream fp) {
        if(p>0) {
            try {
                fp.write(String.format("(%g,%g,%g)", pts[0] * 0.5, pts[1] * 0.5, pts[2] * 0.5).getBytes());
                for (int i = 1; i < p; i++)
                    fp.write(String.format(" (%g,%g,%g)", pts[4 * i + 0] * 0.5, pts[4 * i + 1] * 0.5, pts[4 * i + 2] * 0.5).getBytes());
            }
            catch (IOException ex) {
                throw new RuntimeException(ex);
            }
        }
    }

    /** Returns a vector of the vertex vectors in the global coordinate system.
     * \param[out] v the vector to store the results in.
     * \param[in] (x,y,z) the position vector of the particle in the global
     *                    coordinate system. */
    public void vertices(double x,double y,double z, DoubleArrayList v) {
        v.clear();
        for(int i=0;i<p;i++) {
            v.add(x+pts[4*i+0]*0.5);
            v.add(y+pts[4*i+1]*0.5);
            v.add(z+pts[4*i+2]*0.5);
        }
    }

    /** Outputs the vertex vectors using the global coordinate system.
     * \param[out] fp the file handle to write to.
     * \param[in] (x,y,z) the position vector of the particle in the global
     *                    coordinate system. */
    public void output_vertices(double x, double y, double z) {
        output_vertices(x, y, z, System.out);
    }
    public void output_vertices(double x,double y,double z, OutputStream fp) {
        if(p>0) {
            try {
                fp.write(String.format("(%g,%g,%g)", x + pts[0] * 0.5, y + pts[1] * 0.5, z + pts[2] * 0.5).getBytes());
                for (int i = 1; i < p; i++)
                    fp.write(String.format(" (%g,%g,%g)", x + pts[4 * i + 0] * 0.5, y + pts[4 * i + 1] * 0.5, z + pts[4 * i + 2] * 0.5).getBytes());
            }
            catch (IOException ex) {
                throw new RuntimeException(ex);
            }
        }
    }

    /** Calculates the solid angles of each face of the Voronoi cell and prints
     * the results to an output stream.
     * \param[out] v the vector to store the results in. */
    public void solid_angles(DoubleArrayList v) {
        double solid_angle;
        double[][] normalized = new double[p][3];
        v.clear();
        int i,j,k,l,m,n;
        for (i=0;i<p;i++) {
            normalize_vector(pts, i, normalized[i]);
        }

        for(i=1;i<p;i++) for(j=0;j<nu[i];j++) {
            k=ed_.get(i,j);
            if(k>=0) {
                solid_angle=0;
                ed_.set(i,j,-1-k);
                l=cycle_up(ed_.get(i,nu[i]+j),k);
                m=ed_.get(k,l);
                ed_.set(k,l,-1-m);
                while(m!=i) {
                    n=cycle_up(ed_.get(k,nu[k]+l),m);
                    solid_angle+=calculate_solid_angle(normalized[i],normalized[k],normalized[m]);
                    k=m;l=n;
                    m=ed_.get(k,l);
                    ed_.set(k,l,-1-m);
                }
                v.add(solid_angle);
            }
        }
        reset_edges();
    }

    /** Calculates the areas of each face of the Voronoi cell and prints the
     * results to an output stream.
     * \param[out] v the vector to store the results in. */
    public void face_areas(DoubleArrayList v) {
        double area;
        v.clear();
        int i,j,k,l,m,n;
        double ux,uy,uz,vx,vy,vz,wx,wy,wz;
        for(i=1;i<p;i++) for(j=0;j<nu[i];j++) {
            k = ed_.get(i,j);
            if(k>=0) {
                area=0;
                ed_.set(i,j,-1-k);;
                l = cycle_up(ed_.get(i,nu[i]+j),k);
                m=ed_.get(k,l);ed_.set(k,l,-1-m);
                while(m!=i) {
                    n=cycle_up(ed_.get(k,nu[k]+l),m);
                    ux=pts[4*k]-pts[4*i];
                    uy=pts[4*k+1]-pts[4*i+1];
                    uz=pts[4*k+2]-pts[4*i+2];
                    vx=pts[4*m]-pts[4*i];
                    vy=pts[4*m+1]-pts[4*i+1];
                    vz=pts[4*m+2]-pts[4*i+2];
                    wx=uy*vz-uz*vy;
                    wy=uz*vx-ux*vz;
                    wz=ux*vy-uy*vx;
                    area+=Math.sqrt(wx*wx+wy*wy+wz*wz);
                    k=m;l=n;
                    m=ed_.get(k,l);ed_.set(k,l,-1-m);
                }
                v.add(0.125*area);
            }
        }
        reset_edges();
    }

    /** Calculates the contributions to the Minkowski functionals for this Voronoi cell.
     * \param[in] r the radius to consider.
     * \param[out] ar the area functional.
     * \param[out] vo the volume functional. */
    public void minkowski(double r, double[] ar, double[] vo) {
        int i,j,k,l,m,n;
        ar[0]=vo[0]=0;
        r*=2;
        for(i=1;i<p;i++) for(j=0;j<nu[i];j++) {
            k=ed_.get(i,j);
            if(k>=0) {
                ed_.set(i,j,-1-k);
                l=cycle_up(ed_.get(i,nu[i]+j),k);
                m=ed_.get(k,l);
                ed_.set(k,l,-1-m);
                while(m!=i) {
                    n=cycle_up(ed_.get(k,nu[k]+l),m);
                    minkowski_contrib(i,k,m,r,ar,vo);
                    k=m;l=n;
                    m=ed_.get(k,l);
                    ed_.set(k,l,-1-m);
                }
            }
        }
        vo[0]*=0.125;
        ar[0]*=0.25;
        reset_edges();
    }
    public void output_solid_angles() {
        output_solid_angles(System.out);
    }
    public void output_solid_angles(OutputStream fp) {
        DoubleArrayList v = new DoubleArrayList();
        solid_angles(v);
        Common.voro_print_vector(v, fp);
    }
    /** Outputs the areas of the faces.
     * \param[in] fp the file handle to write to. */
    public void output_face_areas() {
        output_face_areas(System.out);
    }
    public void output_face_areas(OutputStream out) {
        DoubleArrayList v = new DoubleArrayList();
        face_areas(v);
        voro_print_vector(v,out);
    }

    /** Outputs a list of the number of edges in each face.
     * \param[out] v the vector to store the results in. */
    public void face_orders(IntArrayList v) {
        v.clear();
        for(int i=1;i<p;i++)
            for(int j=0;j<nu[i];j++) {
                int k = ed_.get(i,j);
                if(k>=0) {
                    int q=1;
                    ed_.set(i,j,-1-k);;
                    int l = cycle_up(ed_.get(i,nu[i]+j),k);
                    do {
                        q++;
                        int m=ed_.get(k,l);
                        ed_.set(k,l,-1-m);
                        l=cycle_up(ed_.get(k,nu[k]+l),m);
                        k=m;
                    } while (k!=i);
                    v.add(q);;
                }
            }
        reset_edges();
    }
    /** Outputs a list of the number of sides of each face.
     * \param[in] fp the file handle to write to. */
    public void output_face_orders() {
        output_face_orders(System.out);
    }
    public void output_face_orders(OutputStream out) {
        IntArrayList v = new IntArrayList();
        face_orders(v);
        voro_print_vector(v,out);
    }

    /** Computes the number of edges that each face has and outputs a frequency
     * table of the results.
     * \param[out] v the vector to store the results in. */
    public void face_freq_table(IntArrayList v) {
        v.clear();
        for(int i=1;i<p;i++)
            for(int j=0;j<nu[i];j++) {
                int k = ed_.get(i,j);
                if(k>=0) {
                    int q=1;
                    ed_.set(i,j,-1-k);;
                    int l = cycle_up(ed_.get(i,nu[i]+j),k);
                    do {
                        q++;
                        int m=ed_.get(k,l);
                        ed_.set(k,l,-1-m);
                        l=cycle_up(ed_.get(k,nu[k]+l),m);
                        k=m;
                    } while (k!=i);
                    while (q >= v.size()) v.add(0);
                    v.set(q, v.getInt(q)+1);
                }
            }
        reset_edges();
    }
    /** Outputs a */
    public void output_face_freq_table() {
        output_face_freq_table(System.out);
    }
    public  void output_face_freq_table(OutputStream out) {
        IntArrayList v = new IntArrayList();
        face_freq_table(v);
        voro_print_vector(v,out);
    }

    /** For each face, this routine outputs a bracketed sequence of numbers
     * containing a list of all the vertices that make up that face.
     * \param[out] v the vector to store the results in. */
    public void face_vertices(IntArrayList v) {
        int vp=0;
        v.clear();
        for(int i=1;i<p;i++)
            for(int j=0;j<nu[i];j++) {
                int k = ed_.get(i,j);
                if(k>=0) {
                    v.add(0);
                    v.add(i);
                    ed_.set(i,j,-1-k);;
                    int l = cycle_up(ed_.get(i,nu[i]+j),k);
                    do {
                        v.add(k);
                        int m=ed_.get(k,l);
                        ed_.set(k,l,-1-m);
                        l=cycle_up(ed_.get(k,nu[k]+l),m);
                        k=m;
                    } while (k!=i);
                    int vn=v.size();
                    v.set(vp, vn-vp-1);
                    vp=vn;
                }
            }
        reset_edges();
    }
    public void output_face_vertices() {
        output_face_vertices(System.out);
    }
    /** Outputs the */
    public void output_face_vertices(OutputStream out) {
        IntArrayList v = new IntArrayList();
        face_vertices(v);
        voro_print_face_vertices(v,out);
    }

    /** This routine returns the perimeters of each face.
     * \param[out] v the vector to store the results in. */
    public void face_perimeters(DoubleArrayList v) {
        v.clear();
        int i,j,k,l,m;
        double dx,dy,dz,perim;
        for(i=1;i<p;i++) for(j=0;j<nu[i];j++) {
            k = ed_.get(i,j);
            if(k>=0) {
                dx=pts[4*k]-pts[4*i];
                dy=pts[4*k+1]-pts[4*i+1];
                dz=pts[4*k+2]-pts[4*i+2];
                perim=Math.sqrt(dx*dx+dy*dy+dz*dz);
                ed_.set(i,j,-1-k);;
                l = cycle_up(ed_.get(i,nu[i]+j),k);
                do {
                    m=ed_.get(k,l);
                    dx=pts[4*m]-pts[4*k];
                    dy=pts[4*m+1]-pts[4*k+1];
                    dz=pts[4*m+2]-pts[4*k+2];
                    perim+=Math.sqrt(dx*dx+dy*dy+dz*dz);
                    ed_.set(k,l,-1-m);
                    l=cycle_up(ed_.get(k,nu[k]+l),m);
                    k=m;
                } while (k!=i);
                v.add(0.5*perim);
            }
        }
        reset_edges();
    }
    /** Outputs a list of the perimeters of each face.
     * \param[in] fp the file handle to write to. */
    public void output_face_perimeters() {
        output_face_perimeters(System.out);
    }
    public void output_face_perimeters(OutputStream out) {
        DoubleArrayList v = new DoubleArrayList();
        face_perimeters(v);
        voro_print_vector(v,out);
    }

    /** This routine calculates the unit normal vectors for every face.
     * \param[out] v the vector to store the results in. */
    public void normals(DoubleArrayList v) {
        int i,j,k;
        v.clear();
        for(i=1;i<p;i++) for(j=0;j<nu[i];j++) {
            k = ed_.get(i,j);
            if(k>=0) normals_search(v,i,j,k);
        }
        reset_edges();
    }
    /** Outputs a list of the perimeters of each face.
     * \param[in] fp the file handle to write to. */
    public void output_normals() {
        output_normals(System.out);
    }
    public void output_normals(OutputStream out) {
        DoubleArrayList v = new DoubleArrayList();
        normals(v);
        voro_print_positions(v,out);
    }
    /** Outputs a custom string of information about the Voronoi
     * cell to a file. It assumes the cell is at (0,0,0) and has a
     * the default_radius associated with it.
     * \param[in] format the custom format string to use.
     * \param[in] fp the file handle to write to. */
    public void output_custom(String format) {
        output_custom(format, System.out);
    }
    public void output_custom(String format, OutputStream out) {output_custom(format,0,0,0,0,Config.default_radius,out);}

    /** Outputs a custom string of information about the Voronoi cell. The string
     * of information follows a similar style as the C printf command, and detailed
     * information about its format is available at
     * http://math.lbl.gov/voro++/doc/custom.html.
     * \param[in] format the custom string to print.
     * \param[in] i the ID of the particle associated with this Voronoi cell.
     * \param[in] (x,y,z) the position of the particle associated with this Voronoi
     *                    cell.
     * \param[in] r a radius associated with the particle.
     * \param[in] fp the file handle to write to. */
    public void output_custom(String format, int i, double x, double y, double z, double r) {
        output_custom(format, i, x, y, z, r, System.out);
    }
    public void output_custom(String format,int i,double x,double y,double z,double r,OutputStream fp) {
        IntArrayList vi = new IntArrayList();
        DoubleArrayList vd = new DoubleArrayList();
        try {
            for (int j=0; j<format.length(); j++) {
                char fmp = format.charAt(j);
                if(fmp=='%' && j<format.length()-1) {
                    j++;
                    fmp = format.charAt(j);
                    switch(fmp) {

                        // Particle-related output
                        case 'i': fp.write(String.format("%d",i).getBytes());break;
                        case 'x': fp.write(String.format("%g",x).getBytes());break;
                        case 'y': fp.write(String.format("%g",y).getBytes());break;
                        case 'z': fp.write(String.format("%g",z).getBytes());break;
                        case 'q': fp.write(String.format("%g %g %g",x,y,z).getBytes());break;
                        case 'r': fp.write(String.format("%g",r).getBytes());break;

                        // Vertex-related output
                        case 'w': fp.write(String.format("%d",p).getBytes());break;
                        case 'p': output_vertices(fp);break;
                        case 'P': output_vertices(x,y,z,fp);break;
                        case 'o': output_vertex_orders(fp);break;
                        case 'm': fp.write(String.format("%g",0.25*max_radius_squared()).getBytes());break;

                        // Edge-related output
                        case 'g': fp.write(String.format("%d",number_of_edges()).getBytes());break;
                        case 'E': fp.write(String.format("%g",total_edge_distance()).getBytes());break;
                        case 'e': face_perimeters(vd);
                            voro_print_vector(vd,fp);break;

                        // Face-related output
                        case 's': fp.write(String.format("%d",number_of_faces()).getBytes());break;
                        case 'F': fp.write(String.format("%g",surface_area()).getBytes());break;
                        case 'A': {
                            face_freq_table(vi);
                            voro_print_vector(vi,fp);
                        } break;
                        case 'a': face_orders(vi);
                            voro_print_vector(vi,fp);break;
                        case 'f': face_areas(vd);
                            voro_print_vector(vd,fp);break;
                        case 't': {
                            face_vertices(vi);
                            voro_print_face_vertices(vi,fp);
                        } break;
                        case 'l': normals(vd);
                            voro_print_positions(vd,fp);
                            break;
                        case 'n': neighbors(vi);
                            voro_print_vector(vi,fp);
                            break;

                        // Volume-related output
                        case 'v': fp.write(String.format("%g",volume()).getBytes());break;
                        case 'c':
                        case 'C':
                            double[] cxyz = centroid();
                            fp.write(String.format("%g %g %g",cxyz[0],cxyz[1],cxyz[2]).getBytes());
                            break;

                        // The percent sign is not part of a
                        // control sequence
                        default: fp.write(("%"+fmp).getBytes());
                    }
                } else fp.write((""+fmp).getBytes());
            }
            fp.write("\n".getBytes());
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }

    /** Cuts the Voronoi cell by a particle whose center is at a separation of
     * (x,y,z) from the cell center. The value of rsq should be initially set to
     * \f$x^2+y^2+z^2\f$.
     * \param[in] vc a reference to the specialized version of the calling class.
     * \param[in] (x,y,z) the normal vector to the plane.
     * \param[in] rsq the distance along this vector of the plane.
     * \param[in] p_id the plane ID (for neighbor tracking only).
     * \return False if the plane cut deleted the cell entirely, true otherwise. */
    public boolean nplane(double x,double y,double z,double rsq,int p_id) {
        int i,j,cp,qp;
        int[] lp = new int[]{up};
        int[] us=new int[1],ls=new int[1];
        int[] uw = new int[1], lw = new int[1];
        double[] u = new double[1];
        double[] l=new double[1];
        stackp = 0;
        up = 0;

        // Initialize the safe testing routine
        px=x;py=y;pz=z;prsq=rsq;
        maskc+=4;
        if (maskc<4) reset_mask();

        uw[0]=m_test(up,u);
        if(uw[0]==2) {
            if(!search_downward(lw,lp,ls,us,l,u)) return false;
            if(lw[0]==1) {up=lp[0];lp[0]=-1;}
        } else if(uw[0]==0) {
            if(!search_upward(uw,lp,ls,us,l,u)) return true;
            if(uw[0]==1) lp[0]=-1;
        } else {
            lp[0]=-1;
        }

        stackp = stackp2 = stackp3 = 0;

        // Store initial number of vertices
        int op=p;

        if(create_facet(lp[0],ls[0],l[0],us[0],u[0],p_id)) return false;
        int k=0;int xtra=0;
        while(k<stackp3) {
            lp[0]=xse[k];
            k++;
            uw[0]=m_test(lp[0],l);
            for(ls[0]=0;ls[0]<nu[lp[0]];ls[0]++) {
                up=ed_.get(lp[0],ls[0]);

                // Skip if this is a new vertex
                uw[0]=m_test(up,u);
                if(up>=op) continue;

                if(uw[0]==0) {
                    if(u[0]>-big_tol&&ed_.get(up,nu[up]<<1)!=-1) {
                        ed_.set(up,nu[up]<<1,-1);
                        if(stackp3==xse.length) add_memory_xse();
                        xse[stackp3] = up;
                        stackp3++;
                    }
                } else if(uw[0]==1) {

                    // This is a possible facet starting
                    // from a vertex on the cutting plane
                    if(create_facet(-1,0,0,0,u[0],p_id)) return false;
                } else {

                    // This is a new facet
                    us[0]=ed_.get(lp[0], nu[lp[0]]+ls[0]);
                    m_test(lp[0],l);
                    if(create_facet(lp[0],ls[0],l[0],us[0],u[0],p_id)) return false;
                }
            }
            xtra++;
        }

        // Reset back pointers on extra search stack
        for(int dsp=0;dsp<stackp3;dsp++) {
            j=xse[dsp];
            ed_.set(j,nu[j]<<1, j);
        }

        // Delete points: first, remove any duplicates
        for (int dsp=0; dsp<stackp; ) {
            j=ds[dsp];
            if(ed_.get(j,nu[j])!=-1) {
                ed_.set(j,nu[j], -1);
                dsp++;
            } else {
                stackp--;
                ds[dsp] = ds[stackp];
            }
        }

        // Add the points in the auxiliary delete stack,
        // and reset their back pointers
        for(int dsp=0;dsp<stackp2;dsp++) {
            j=ds2[dsp];
            ed_.set(j,nu[j]<<1, j);
            if(ed_.get(j,nu[j])!=-1) {
                ed_.set(j,nu[j],-1);
                if(stackp==ds.length) add_memory_ds();
                ds[stackp] = j;
                stackp++;
            }
        }

        // Scan connections and add in extras
        for(int dsp=0;dsp<stackp;dsp++) {
            cp=ds[dsp];
            for(int edp=0;edp<nu[cp];edp++) {
                qp=ed_.get(cp, edp);
                if(qp!=-1&&ed_.get(qp,nu[qp])!=-1) {
                    if(stackp==ds.length) {
                        int dis=stackp-dsp;
                        add_memory_ds();
                        dsp=dis;
                    }
                    ds[stackp] = qp;
                    stackp++;
                    ed_.set(qp,nu[qp],-1);
                }
            }
        }
        up=0;

        // Delete them from the array structure
        while(stackp>0) {
            --p;
            while(ed_.get(p,nu[p])==-1) {
                j=nu[p];
                int edp=0; // ed[p]
                mec[j]--;
                int edd=((j<<1)+1)*mec[j]; // mep[j]
                while(edp<(j<<1)+1) {
                    ed_.set(p,edp,mep[j][edd]);
                    edp++;
                    edd++;
                }
                n_set_aux2_copy(p,j);
                n_copy_pointer(ed_.get(p,j<<1),p);
                ed_.setIndex(ed_.get(p,j<<1), ed_.getStorage(p), ed_.getOffset(p));
                --p;
            }
            stackp--;
            up=ds[stackp];
            if(up<p) {

                // Vertex management
                pts[(up<<2)]=pts[(p<<2)];
                pts[(up<<2)+1]=pts[(p<<2)+1];
                pts[(up<<2)+2]=pts[(p<<2)+2];

                // Memory management
                j=nu[up];
                int edp = 0; // ed[up]
                mec[j]--;
                int edd = ((j<<1)+1)*mec[j]; // mep[j]
                while(edp<(j<<1)+1) {
                    ed_.set(up,edp,mep[j][edd]);
                    edp++;
                    edd++;
                }
                n_set_aux2_copy(up,j);
                n_copy_pointer(ed_.get(up,j<<1),up);
                n_copy_pointer(up,p);

                ed_.setIndex(ed_.get(up,j<<1), ed_.getStorage(up), ed_.getOffset(up));

                // Edge management
                ed_.setIndex(up, ed_.getStorage(p), ed_.getOffset(p));
                nu[up]=nu[p];
                for(i=0;i<nu[up];i++) {
                    ed_.set(ed_.get(up,i), ed_.get(up,nu[up]+i), up);
                }
                ed_.set(up, nu[up]<<1, up);
            } else up=p++;
        }

        // Check for any vertices of zero order
        if(mec[0]>0) voro_fatal_error("Zero order vertex formed",Config.Voropp.INTERNAL_ERROR);
        // Collapse any order 2 vertices and exit
        return collapse_order2();
    }

    /** This routine tests to see whether the cell intersects a plane by starting
     * from the guess point up. If up intersects, then it immediately returns true.
     * Otherwise, it calls the plane_intersects_track() routine.
     * \param[in] (x,y,z) the normal vector to the plane.
     * \param[in] rsq the distance along this vector of the plane.
     * \return False if the plane does not intersect the plane, true if it does. */
    public boolean plane_intersects(double x,double y,double z,double rsq) {
        double g=x*pts[4*up]+y*pts[4*up+1]+z*pts[4*up+2];
        if(g<rsq) return plane_intersects_track(x,y,z,rsq,g);
        return true;
    }

    /** This routine tests to see if a cell intersects a plane. It first tests a
     * random sample of approximately sqrt(p)/4 points. If any of those are
     * intersect, then it immediately returns true. Otherwise, it takes the closest
     * point and passes that to plane_intersect_track() routine.
     * \param[in] (x,y,z) the normal vector to the plane.
     * \param[in] rsq the distance along this vector of the plane.
     * \return False if the plane does not intersect the plane, true if it does. */
    public boolean plane_intersects_guess(double x,double y,double z,double rsq) {
        up=0;
        double g=x*pts[4*up]+y*pts[4*up+1]+z*pts[4*up+2];
        if(g<rsq) {
            int ca=1,cc=p>>3,mp=1;
            double m;
            while(ca<cc) {
                m=x*pts[4*mp]+y*pts[4*mp+1]+z*pts[4*mp+2];
                if(m>g) {
                    if(m>rsq) return true;
                    g=m;up=mp;
                }
                ca+=mp++;
            }
            return plane_intersects_track(x,y,z,rsq,g);
        }
        return true;
    }

    /** Constructs the relational table if the edges have been specified. */
    public void construct_relations() {
        int i,j,k,l;
        for(i=0;i<p;i++) for(j=0;j<nu[i];j++) {
            k = ed_.get(i,j);
            l=0;
            while(ed_.get(k,l)!=i) {
                l++;
                if(l==nu[k]) voro_fatal_error("Relation table construction failed",Config.Voropp.INTERNAL_ERROR);
            }
            ed_.set(i,nu[i]+j,l);
        }
    }

    /** Checks that the relational table of the Voronoi cell is accurate, and
     * prints out any errors. This algorithm is O(p), so running it every time the
     * plane routine is called will result in a significant slowdown. */
    public void check_relations() {
        int i,j;
        for(i=0;i<p;i++) for(j=0;j<nu[i];j++) if(ed_.get(ed_.get(i,j),ed_.get(i,nu[i]+j))!=i)
            System.out.printf("Relational error at point %d, edge %d.\n",i,j);
    }

    /** This routine checks for any two vertices that are connected by more than
     * one edge. The plane algorithm is designed so that this should not happen, so
     * any occurrences are most likely errors. Note that the routine is O(p), so
     * running it every time the plane routine is called will result in a
     * significant slowdown. */
    public void check_duplicates() {
        int i,j,k;
        for(i=0;i<p;i++) for(j=1;j<nu[i];j++) for(k=0;k<j;k++) if(ed_.get(i,j)==ed_.get(i,k))
            System.out.printf("Duplicate edges: (%d,%d) and (%d,%d) [%d]\n",i,j,i,k,ed_.get(i,j));
    }

    /** Prints the vertices, their edges, the relation table, and also notifies if
     * any memory errors are visible. */
    public void print_edges() {
        int j;
        for(int i=0;i<p;i++) {
            System.out.printf("%d %d  ",i,nu[i]);
            for(j=0;j<nu[i];j++) System.out.printf(" %d",ed_.get(i,j));
            System.out.print("  ");
            while(j<(nu[i]<<1)) System.out.printf(" %d",ed_.get(i,j));
            System.out.printf("   %d",ed_.get(i,j));
            print_edges_neighbors(i);
            System.out.printf("  %g %g %g",pts[4*i+0],pts[4*i+1],pts[4*i+2]);
            if (ed_.storage[i] != mep[nu[i]]) throw new RuntimeException("unexpected storage");
            if (ed_.getOffset(i) > mec[nu[i]]*((nu[i]<<1)+1)) System.out.println(" Memory error");
            else System.out.println();
        }
    }
    /** Returns a list of IDs of neighboring particles
     * corresponding to each face.
     * \param[out] v a reference to a vector in which to return the
     *               results. If no neighbor information is
     *               available, a blank vector is returned. */
    public void neighbors(IntArrayList v) {v.clear();}
    /** This is a virtual function that is overridden by a routine
     * to print a list of IDs of neighboring particles
     * corresponding to each face. By default, when no neighbor
     * information is available, the routine does nothing.
     * \param[in] fp the file handle to write to. */
    public void output_neighbors() {
        output_neighbors(System.out);
    }
    public void output_neighbors(OutputStream fp) {}
    /** This a virtual function that is overridden by a routine to
     * print the neighboring particle IDs for a given vertex. By
     * default, when no neighbor information is available, the
     * routine does nothing.
     * \param[in] i the vertex to consider. */
    public void print_edges_neighbors(int i) {}
    /** This is a simple inline function for picking out the index
     * of the next edge counterclockwise at the current vertex.
     * \param[in] a the index of an edge of the current vertex.
     * \param[in] p the number of the vertex.
     * \return 0 if a=nu[p]-1, or a+1 otherwise. */
    public int cycle_up(int a,int p) {return a==nu[p]-1?0:a+1;}
    /** This is a simple inline function for picking out the index
     * of the next edge clockwise from the current vertex.
     * \param[in] a the index of an edge of the current vertex.
     * \param[in] p the number of the vertex.
     * \return nu[p]-1 if a=0, or a-1 otherwise. */
    public int cycle_down(int a,int p) {return a==0?nu[p]-1:a-1;}

    /** This a one dimensional array that holds the current sizes
     * of the memory allocations for them mep array.*/
    protected int[] mem;
    /** This is a one dimensional array that holds the current
     * number of vertices of order p that are stored in the mep[p]
     * array. */
    protected int[] mec;
    /** This is a two dimensional array for holding the information
     * about the edges of the Voronoi cell. mep[p] is a
     * one-dimensional array for holding the edge information about
     * all vertices of order p, with each vertex holding 2*p+1
     * integers of information. The total number of vertices held
     * on mep[p] is stored in mem[p]. If the space runs out, the
     * code allocates more using the add_memory() routine. */
    protected int[][] mep;
    protected void reset_edges() {
        int i,j;
        for(i=0;i<p;i++) for(j=0;j<nu[i];j++) {
            if(ed_.get(i,j)>=0) voro_fatal_error("Edge reset routine found a previously untested edge",Config.Voropp.INTERNAL_ERROR);
            ed_.set(i,j,-1-ed_.get(i,j));
        }
    }
    /** Ensures that enough memory is allocated prior to carrying out a copy.
     * \param[in] vc a reference to the specialized version of the calling class.
     * \param[in] vb a pointered to the class to be copied. */
    protected void check_memory_for_copy(VoronoiCellBase vb) {
        while(current_vertex_order < vb.current_vertex_order) add_memory_vorder();
        for(int i=0;i<current_vertex_order;i++) while(mem[i]<vb.mec[i]) add_memory(i);
        while(current_vertices < vb.p) add_memory_vertices();
    }

    /** Copies the vertex and edge information from another class. The routine
     * assumes that enough memory is available for the copy.
     * \param[in] vb a pointer to the class to copy. */
    protected void copy(VoronoiCellBase vb) {
        int i,j;
        p=vb.p;up=0;
        for(i=0;i<current_vertex_order;i++) {
            mec[i]=vb.mec[i];
            for(j=0;j<mec[i]*(2*i+1);j++) mep[i][j]=vb.mep[i][j];
            for(j=0;j<mec[i]*(2*i+1);j+=2*i+1) ed_.setIndex(mep[i][j+2*i], mep[i], j);
        }
        for(i=0;i<p;i++) nu[i]=vb.nu[i];
        for(i=0;i<(p<<2);i++) pts[i]=vb.pts[i];
    }

    /** This is the delete stack, used to store the vertices which
     * are going to be deleted during the plane cutting procedure.
     */
    private int[] ds;
    private int stackp;
    /** This is the auxiliary delete stack, which has size set by
     * current_delete2_size. */
    private int[] ds2;
    private int stackp2;
    /** This is the exrta search stack. */
    private int[] xse;
    private int stackp3;
    private int maskc;
    /** The x coordinate of the normal vector to the test plane. */
    private double px;
    /** The y coordinate of the normal vector to the test plane. */
    private double py;
    /** The z coordinate of the normal vector to the test plane. */
    private double pz;
    /** The magnitude of the normal vector to the test plane. */
    private double prsq;
    protected abstract void n_allocate(int i, int m);
    protected abstract void n_allocate_aux1(int i);
    protected abstract void n_set_to_aux1_offset(int k, int m);
    protected abstract void n_copy_to_aux1(int i, int m);
    protected abstract void n_switch_to_aux1(int i);
    protected abstract void n_add_memory_vertices(int i);
    protected abstract void n_add_memory_vorder(int i);
    protected abstract void n_set_pointer(int p,int n);
    protected abstract void n_copy(int a,int b,int c,int d);
    protected abstract void n_set(int a,int b,int c);
    protected abstract void n_set_aux1(int k);
    protected abstract void n_copy_aux1(int a,int b);
    protected abstract void n_set_aux2_copy(int a,int b);
    protected abstract void n_copy_pointer(int a,int b);
    protected abstract void n_set_to_aux1(int j);
    protected abstract void n_copy_aux1_shift(int a,int b);
    protected abstract void n_set_to_aux2(int j);

    /** Increases the memory storage for a particular vertex order, by increasing
     * the size of the of the corresponding mep array. If the arrays already exist,
     * their size is doubled; if they don't exist, then new ones of size
     * init_n_vertices are allocated. The routine also ensures that the pointers in
     * the ed array are updated, by making use of the back pointers. For the cases
     * where the back pointer has been temporarily overwritten in the marginal
     * vertex code, the auxiliary delete stack is scanned to find out how to update
     * the ed value. If the template has been instantiated with the neighbor
     * tracking turned on, then the routine also reallocates the corresponding mne
     * array.
     * \param[in] i the order of the vertex memory to be increased. */
    private void add_memory(int i) {
        int s=(i<<1)+1;
        if(mem[i]==0) {
            n_allocate(i,Config.init_n_vertices);
            mep[i]=new int[Config.init_n_vertices*s];
            mem[i]=Config.init_n_vertices;
            if (Config.VOROPP_VERBOSE >=2) {
                System.err.printf("Order %d vertex memory created\n", i);
            }
        } else {
            int j=0,k;
            mem[i]<<=1;
            if(mem[i]>Config.max_n_vertices) voro_fatal_error("Point memory allocation exceeded absolute maximum", Config.Voropp.MEMORY_ERROR);
            if (Config.VOROPP_VERBOSE >= 2) {
                System.err.printf("Order %d vertex memory scaled up to %d\n", i, mem[i]);
            }
            int[] l=new int[s*mem[i]];
            boolean[] migrated = new boolean[ed_.offsets.length];
            int m=0;
            n_allocate_aux1(i);
            while(j<s*mec[i]) {
                k=mep[i][j+(i<<1)];
                if(k>=0) {
                    ed_.setIndex(k, l, j);
                    migrated[k] = true;
                    n_set_to_aux1_offset(k,m);
                } else {
                    int dsp;
                    for(dsp=0;dsp<stackp2;dsp++) {
                        if(!migrated[ds2[dsp]]) {
                            ed_.setIndex(ds2[dsp], l, j);
                            n_set_to_aux1_offset(ds2[dsp],m);
                            break;
                        }
                    }
                    if(dsp==stackp2) {
                        for (dsp=0; dsp<stackp3; dsp++) {
                            if (ed_.getStorage(xse[dsp]) == mep[i] && ed_.getOffset(xse[dsp]) == j) {

                                ed_.setIndex(xse[dsp], l, j);
                                n_set_to_aux1_offset(xse[dsp], m);
                                break;
                            }
                        }
                        if (dsp == stackp3) voro_fatal_error("Couldn't relocate dangling pointer",Config.Voropp.INTERNAL_ERROR);
                    }
                    if (Config.VOROPP_VERBOSE >= 3) {
                        System.err.println("Relocated dangling pointer");
                    }
                }
                for(k=0;k<s;k++,j++) l[j]=mep[i][j];
                for(k=0;k<i;k++,m++) n_copy_to_aux1(i,m);
            }
            mep[i]=l;
            n_switch_to_aux1(i);
        }

    }

    /** Doubles the maximum number of vertices allowed, by reallocating the ed, nu,
     * and pts arrays. If the allocation exceeds the absolute maximum set in
     * max_vertices, then the routine exits with a fatal error. If the template has
     * been instantiated with the neighbor tracking turned on, then the routine
     * also reallocates the ne array. */
    private void add_memory_vertices() {
        int i=(current_vertices<<1);
        if(i>Config.max_vertices) voro_fatal_error("Vertex memory allocation exceeded absolute maximum",Config.Voropp.MEMORY_ERROR);
        if (Config.VOROPP_VERBOSE >= 2) {
            System.err.printf("Vertex memory scaled up to %d\n", i);
        }

        n_add_memory_vertices(i);
        nu = Arrays.copyOf(nu, i);
        mask = Arrays.copyOf(mask, i);
        pts = Arrays.copyOf(pts, 4*i);
        current_vertices=i;

    }

    /** Doubles the maximum allowed vertex order, by reallocating mem, mep, and mec
     * arrays. If the allocation exceeds the absolute maximum set in
     * max_vertex_order, then the routine causes a fatal error. If the template has
     * been instantiated with the neighbor tracking turned on, then the routine
     * also reallocates the mne array. */
    private void add_memory_vorder() {
        int i=(current_vertex_order<<1),j;
        if(i>Config.max_vertex_order) voro_fatal_error("Vertex order memory allocation exceeded absolute maximum",Config.Voropp.MEMORY_ERROR);
        if (Config.VOROPP_VERBOSE >= 2) {
            System.err.printf("Vertex order memory scaled up to %d\n", i);
        }
        mem = Arrays.copyOf(mem, i);
        mep = Arrays.copyOf(mep, i);
        mec = Arrays.copyOf(mec, i);
        n_add_memory_vorder(i);
        current_vertex_order=i;

    }

    /** Doubles the size allocation of the main delete stack. If the allocation
     * exceeds the absolute maximum set in max_delete_size, then routine causes a
     * fatal error. */
    private void add_memory_ds() {
        current_delete_size<<=1;
        if(current_delete_size>Config.max_delete_size) voro_fatal_error("Delete stack 1 memory allocation exceeded absolute maximum",Config.Voropp.MEMORY_ERROR);
        if (Config.VOROPP_VERBOSE >= 2) {
            System.err.printf("Delete stack 1 memory scaled up to %d\n", current_delete_size);
        }
        ds = Arrays.copyOf(ds, current_delete_size);
    }

    /** Doubles the size allocation of the auxiliary delete stack. If the
     * allocation exceeds the absolute maximum set in max_delete2_size, then the
     * routine causes a fatal error. */
    private void add_memory_ds2() {
        current_delete2_size<<=1;
        if(current_delete2_size>Config.max_delete2_size) voro_fatal_error("Delete stack 2 memory allocation exceeded absolute maximum",Config.Voropp.MEMORY_ERROR);
        if (Config.VOROPP_VERBOSE >= 2) {
            System.err.printf("Delete stack 2 memory scaled up to %d\n", current_delete2_size);
        }
        ds2 = Arrays.copyOf(ds2, current_delete2_size);
    }

    /** Doubles the size allocation of the auxiliary delete stack. If the
     * allocation exceeds the absolute maximum set in max_delete2_size, then the
     * routine causes a fatal error. */
    private void add_memory_xse() {
        current_xsearch_size<<=1;
        if(current_xsearch_size>Config.max_xsearch_size) voro_fatal_error("Extra search stack memory allocation exceeded absolute maximum",Config.Voropp.MEMORY_ERROR);
        if (Config.VOROPP_VERBOSE >= 2) {
            System.err.printf("Extra search stack memory scaled up to %d\n", current_xsearch_size);
        }
        xse = Arrays.copyOf(xse, current_xsearch_size);
    }

    /** This routine is a fall-back, in case floating point errors caused the usual
     * search routine to fail. In the fall-back routine, we just test every edge to
     * find one straddling the plane. */
    private boolean failsafe_find(int[] lp, int[] ls, int[] us, double[] l, double[] u) {
        throw new MethodNotImplementedException("Bailed out of convex calculation (not supported yet)");
    }

    /** Creates a new facet.
     * \return True if cell deleted, false otherwise. */
    private boolean create_facet(int lp, int ls, double l, int us, double u, int p_id) {
        int i,j,k,qp,qs,iqs,cp,cs,rp,edp,edd;
        int lw,qw;
        boolean new_double_edge=false,double_edge=false;
        double q=0,r;

        // We're about to add the first point of the new facet. In either
        // routine, we have to add a point, so first check there's space for
        // it.
        if(p==current_vertices) add_memory_vertices();

        if(lp==-1) {

            // We want to be strict about reaching the conclusion that the
            // cell is entirely within the cutting plane. It's not enough
            // to find a vertex that has edges which are all inside or on
            // the plane. If the vertex has neighbors that are also on the
            // plane, we should check those too.
            int[] upout = new int[]{up};
            if(!search_for_outside_edge(upout)) {
                up = upout[0];
                return true;
            }
            up = upout[0];

            // The search algorithm found a point which is on the cutting
            // plane. We leave that point in place, and create a new one at
            // the same location.
            pts[(p<<2)]=pts[(up<<2)];
            pts[(p<<2)+1]=pts[(up<<2)+1];
            pts[(p<<2)+2]=pts[(up<<2)+2];

            // Search for a collection of edges of the test vertex which
            // are outside of the cutting space. Begin by testing the
            // zeroth edge.
            i=0;
            lp=ed_.get(up,0);
            double[] lout = new double[]{l};
            lw=m_testx(lp,lout);
            l = lout[0];
            if(lw!=0) {

                // The first edge is either inside the cutting space,
                // or lies within the cutting plane. Test the edges
                // sequentially until we find one that is outside.
                int rw=lw;
                do {
                    i++;

                    // If we reached the last edge with no luck
                    // then all of the vertices are inside
                    // or on the plane, so the cell is completely
                    // deleted
                    if(i==nu[up]) return true;
                    lp=ed_.get(up,i);
                    lout = new double[]{l};
                    lw=m_testx(lp,lout);
                    l = lout[0];
                } while (lw!=0);
                j=i+1;

                // We found an edge outside the cutting space. Keep
                // moving through these edges until we find one that's
                // inside or on the plane.
                while(j<nu[up]) {
                    lp=ed_.get(up,j);
                    lout = new double[]{l};
                    lw=m_testx(lp,lout);
                    l = lout[0];
                    if(lw!=0) break;
                    j++;
                }

                // Compute the number of edges for the new vertex. In
                // general it will be the number of outside edges
                // found, plus two. But we need to recognize the
                // special case when all but one edge is outside, and
                // the remaining one is on the plane. For that case we
                // have to reduce the edge count by one to prevent
                // doubling up.
                if(j==nu[up]&&i==1&&rw==1) {
                    nu[p]=nu[up];
                    double_edge=true;
                } else nu[p]=j-i+2;
                k=1;

                // Add memory for the new vertex if needed, and
                // initialize
                while (nu[p]>=current_vertex_order) add_memory_vorder();
                if(mec[nu[p]]==mem[nu[p]]) add_memory(nu[p]);
                n_set_pointer(p,nu[p]);
                ed_.setIndex(p, mep[nu[p]], ((nu[p]<<1)+1)*mec[nu[p]]);
                mec[nu[p]]++;
                ed_.set(p, nu[p]<<1, p);

                // Copy the edges of the original vertex into the new
                // one. Delete the edges of the original vertex, and
                // update the relational table.
                us=cycle_down(i,up);
                while(i<j) {
                    qp=ed_.get(up,i);
                    qs=ed_.get(up,nu[up]+i);
                    n_copy(p,k,up,i);
                    ed_.set(p,k,qp);
                    ed_.set(p,nu[p]+k,qs);
                    ed_.set(qp,qs,p);
                    ed_.set(qp,nu[qp]+qs,k);
                    ed_.set(up,i,-1);
                    i++;k++;
                }
                qs=i==nu[up]?0:i;
            } else {

                // In this case, the zeroth edge is outside the cutting
                // plane. Begin by searching backwards from the last
                // edge until we find an edge which isn't outside.
                i=nu[up]-1;
                lp=ed_.get(up,i);
                lout = new double[]{l};
                lw=m_testx(lp,lout);
                l = lout[0];
                while(lw==0) {
                    i--;

                    // If i reaches zero, then we have a point in
                    // the plane all of whose edges are outside
                    // the cutting space, so we just exit
                    if(i==0) return false;
                    lp=ed_.get(up,i);
                    lout = new double[]{l};
                    lw=m_testx(lp,lout);
                    l = lout[0];
                }

                // Now search forwards from zero
                j=1;
                qp=ed_.get(up,j);
                double[] qout = new double[]{q};
                qw=m_testx(qp,qout);
                q = qout[0];
                while(qw==0) {
                    j++;
                    qp=ed_.get(up,j);
                    lout = new double[]{l};
                    qw=m_testx(qp,lout);
                    l = lout[0];
                }

                // Compute the number of edges for the new vertex. In
                // general it will be the number of outside edges
                // found, plus two. But we need to recognize the
                // special case when all but one edge is outside, and
                // the remaining one is on the plane. For that case we
                // have to reduce the edge count by one to prevent
                // doubling up.
                if(i==j&&qw==1) {
                    double_edge=true;
                    nu[p]=nu[up];
                } else {
                    nu[p]=nu[up]-i+j+1;
                }

                // Add memory to store the vertex if it doesn't exist
                // already
                k=1;
                while(nu[p]>=current_vertex_order) add_memory_vorder();
                if(mec[nu[p]]==mem[nu[p]]) add_memory(nu[p]);

                // Copy the edges of the original vertex into the new
                // one. Delete the edges of the original vertex, and
                // update the relational table.
                n_set_pointer(p,nu[p]);
                ed_.setIndex(p, mep[nu[p]], ((nu[p]<<1)+1)*mec[nu[p]]);
                mec[nu[p]]++;
                ed_.set(p, nu[p]<<1, p);
                us=i++;
                while(i<nu[up]) {
                    qp=ed_.get(up,i);
                    qs=ed_.get(up,nu[up]+i);
                    n_copy(p,k,up,i);
                    ed_.set(p,k,qp);
                    ed_.set(p,nu[p]+k,qs);
                    ed_.set(qp,qs,p);
                    ed_.set(qp,nu[qp]+qs,k);
                    ed_.set(up,i,-1);
                    i++;k++;
                }
                i=0;
                while(i<j) {
                    qp=ed_.get(up,i);
                    qs=ed_.get(up,nu[up]+i);
                    n_copy(p,k,up,i);
                    ed_.set(p,k,qp);
                    ed_.set(p,nu[p]+k,qs);
                    ed_.set(qp,qs,p);
                    ed_.set(qp,nu[qp]+qs,k);
                    ed_.set(up,i,-1);
                    i++;
                    k++;
                }
                qs=j;
            }
            if(!double_edge) {
                n_copy(p,k,up,qs);
                n_set(p,0,p_id);
            } else n_copy(p,0,up,qs);

            // Add this point to the auxiliary delete stack
            if(stackp2==ds2.length) add_memory_ds2();
            ds2[stackp2] = up;
            stackp2++;

            // Look at the edges on either side of the group that was
            // detected. We're going to commence facet computation by
            // moving along one of them. We are going to end up coming back
            // along the other one.
            cs=k;
            qp=up;q=u;
            i=ed_.get(up,us);
            us=ed_.get(up,nu[up]+us);
            up=i;
            ed_.set(qp,nu[qp]<<1,-p);

        } else {

            // The search algorithm found an intersected edge between the
            // points lp and up. Create a new vertex between them which
            // lies on the cutting plane. Since u and l differ by at least
            // the tolerance, this division should never screw up.
            if(stackp==ds.length) add_memory_ds();
            ds[stackp] = up;
            stackp++;
            r=u/(u-l);l=1-r;
            pts[p<<2]=pts[lp<<2]*r+pts[up<<2]*l;
            pts[(p<<2)+1]=pts[(lp<<2)+1]*r+pts[(up<<2)+1]*l;
            pts[(p<<2)+2]=pts[(lp<<2)+2]*r+pts[(up<<2)+2]*l;

            // This point will always have three edges. Connect one of them
            // to lp.
            nu[p]=3;
            if(mec[3]==mem[3]) add_memory(3);
            n_set_pointer(p,3);
            n_set(p,0,p_id);
            n_copy(p,1,up,us);
            n_copy(p,2,lp,ls);
            ed_.setIndex(p, mep[3], 7*mec[3]);
            mec[3]++;
            ed_.set(p, 6, p);
            ed_.set(up,us, -1);
            ed_.set(lp,ls,p);
            ed_.set(lp,nu[lp]+ls,1);
            ed_.set(p, 1, lp);
            ed_.set(p, nu[p]+1, ls);
            cs=2;

            // Set the direction to move in
            qs=cycle_up(us,up);
            qp=up;q=u;
        }

        // When the code reaches here, we have initialized the first point, and
        // we have a direction for moving it to construct the rest of the facet
        cp=p;rp=p;p++;
        while(qp!=up||qs!=us) {

            // We're currently tracing round an intersected facet. Keep
            // moving around it until we find a point or edge which
            // intersects the plane.
            lp=ed_.get(qp,qs);
            double[] lout = new double[]{l};
            lw=m_testx(lp,lout);
            l = lout[0];

            if(lw==2) {

                // The point is still in the cutting space. Just add it
                // to the delete stack and keep moving.
                qs=cycle_up(ed_.get(qp,nu[qp]+qs),lp);
                qp=lp;
                q=l;
                if(stackp==ds.length) add_memory_ds();
                ds[stackp] = qp;
                stackp++;

            } else if(lw==0) {

                // The point is outside of the cutting space, so we've
                // found an intersected edge. Introduce a regular point
                // at the point of intersection. Connect it to the
                // point we just tested. Also connect it to the previous
                // new point in the facet we're constructing.
                if(p==current_vertices) add_memory_vertices();
                r=q/(q-l);l=1-r;
                pts[4*p]=pts[4*lp]*r+pts[4*qp]*l;
                pts[4*p+1]=pts[4*lp+1]*r+pts[4*qp+1]*l;
                pts[4*p+2]=pts[4*lp+2]*r+pts[4*qp+2]*l;
                nu[p]=3;
                if(mec[3]==mem[3]) add_memory(3);
                ls=ed_.get(qp,qs+nu[qp]);
                n_set_pointer(p,3);
                n_set(p,0,p_id);
                n_copy(p,1,qp,qs);
                n_copy(p,2,lp,ls);
                ed_.setIndex(p, mep[3], 7*mec[3]);
                mec[3]++;
			    ed_.set(p,0,cp);
                ed_.set(p,1,lp);
                ed_.set(p,3,cs);
                ed_.set(p,4,ls);
                ed_.set(p,6,p);
                ed_.set(lp,ls,p);
                ed_.set(lp,nu[lp]+ls,1);
                ed_.set(cp,cs,p);
                ed_.set(cp,nu[cp]+cs,0);
                ed_.set(qp,qs,-1);
                qs=cycle_up(qs,qp);
                cp=p++;
                cs=2;
            } else {

                // We've found a point which is on the cutting plane.
                // We're going to introduce a new point right here, but
                // first we need to figure out the number of edges it
                // has.
                if(p==current_vertices) add_memory_vertices();

                // If the previous vertex detected a double edge, our
                // new vertex will have one less edge.
                k=double_edge?0:1;
                qs=ed_.get(qp,nu[qp]+qs);
                qp=lp;
                iqs=qs;

                // Start testing the edges of the current point until
                // we find one which isn't outside the cutting space
                do {
                    k++;
                    qs=cycle_up(qs,qp);
                    lp=ed_.get(qp,qs);
                    lout = new double[]{l};
                    lw=m_testx(lp,lout);
                    l=lout[0];
                } while (lw==0);

                // Now we need to find out whether this marginal vertex
                // we are on has been visited before, because if that's
                // the case, we need to add vertices to the existing
                // new vertex, rather than creating a fresh one. We also
                // need to figure out whether we're in a case where we
                // might be creating a duplicate edge.
                j=-ed_.get(qp,nu[qp]<<1);
                if(qp==up&&qs==us) {

                    // If we're heading into the final part of the
                    // new facet, then we never worry about the
                    // duplicate edge calculation.
                    new_double_edge=false;
                    if(j>0) k+=nu[j];
                } else {
                    if(j>0) {

                        // This vertex was visited before, so
                        // count those vertices to the ones we
                        // already have.
                        k+=nu[j];

                        // The only time when we might make a
                        // duplicate edge is if the point we're
                        // going to move to next is also a
                        // marginal point, so test for that
                        // first.
                        if(lw==1) {

                            // Now see whether this marginal point
                            // has been visited before.
                            i=-ed_.get(lp,nu[lp]<<1);
                            if(i>0) {

                                // Now see if the last edge of that other
                                // marginal point actually ends up here.
                                if(ed_.get(i,nu[i]-1)==j) {
                                    new_double_edge=true;
                                    k-=1;
                                } else new_double_edge=false;
                            } else {

                                // That marginal point hasn't been visited
                                // before, so we probably don't have to worry
                                // about duplicate edges, except in the
                                // case when that's the way into the end
                                // of the facet, because that way always creates
                                // an edge.
                                if(j==rp&&lp==up&&ed_.get(qp,nu[qp]+qs)==us) {
                                    new_double_edge=true;
                                    k-=1;
                                } else new_double_edge=false;
                            }
                        } else new_double_edge=false;
                    } else {

                        // The vertex hasn't been visited
                        // before, but let's see if it's
                        // marginal
                        if(lw==1) {

                            // If it is, we need to check
                            // for the case that it's a
                            // small branch, and that we're
                            // heading right back to where
                            // we came from
                            i=-ed_.get(lp,nu[lp]<<1);
                            if(i==cp) {
                                new_double_edge=true;
                                k-=1;
                            } else new_double_edge=false;
                        } else new_double_edge=false;
                    }
                }

                // k now holds the number of edges of the new vertex
                // we are forming. Add memory for it if it doesn't exist
                // already.
                while(k>=current_vertex_order) add_memory_vorder();
                if(mec[k]==mem[k]) add_memory(k);

                // Now create a new vertex with order k, or augment
                // the existing one
                if(j>0) {

                    // If we're augmenting a vertex but we don't
                    // actually need any more edges, just skip this
                    // routine to avoid memory confusion
                    if(nu[j]!=k) {

                        // Allocate memory and copy the edges
                        // of the previous instance into it
                        n_set_aux1(k);
                        edp = ((k<<1)+1)*mec[k]; // mep[k]
                        mec[k]++;
                        i=0;
                        while(i<nu[j]) {
                            n_copy_aux1(j,i);
                            mep[k][edp+i]=ed_.get(j,i);
                            mep[k][edp+k+i]=ed_.get(j,nu[j]+i);
                            i++;
                        }
                        mep[k][edp + k<<1]=j;

                        // Remove the previous instance with
                        // fewer vertices from the memory
                        // structure
                        mec[nu[j]]--;
                        edd = ((nu[j]<<1)+1)*mec[nu[j]]; // mep[nu[j]]
                        if(mep[nu[j]]!=ed_.getStorage(j) || edd != ed_.getOffset(j)) {
                            for(int lll=0;lll<=(nu[j]<<1);lll++) ed_.set(j,lll,mep[nu[j]][edd+lll]);
                            n_set_aux2_copy(j,nu[j]);
                            n_copy_pointer(edd + nu[j]<<1,j); // ??????
                            ed_.setIndex(mep[nu[j]][edd+nu[j]<<1], ed_.getStorage(j), 0);
                        }
                        n_set_to_aux1(j);
                        ed_.setIndex(j, mep[k], edp);
                    } else i=nu[j];
                } else {

                    // Allocate a new vertex of order k
                    n_set_pointer(p,k);
                    ed_.setIndex(p, mep[k], ((k<<1)+1)*mec[k]);
                    mec[k]++;
                    ed_.set(p,k<<1,p);
                    if(stackp2==ds2.length) add_memory_ds2();
                    ds2[stackp2] = qp;
                    stackp2++;
                    pts[4*p]=pts[4*qp];
                    pts[4*p+1]=pts[4*qp+1];
                    pts[4*p+2]=pts[4*qp+2];
                    ed_.set(qp,nu[qp]<<1,-p);
                    j=p++;
                    i=0;
                }
                nu[j]=k;

                // Unless the previous case was a double edge, connect
                // the first available edge of the new vertex to the
                // last one in the facet
                if(!double_edge) {
                    ed_.set(j,i,cp);
                    ed_.set(j,nu[j]+i,cs);
                    n_set(j,i,p_id);
                    ed_.set(cp,cs,j);
                    ed_.set(cp,nu[cp]+cs,i);
                    i++;
                }

                // Copy in the edges of the underlying vertex,
                // and do one less if this was a double edge
                qs=iqs;
                while(i<(new_double_edge?k:k-1)) {
                    qs=cycle_up(qs,qp);
                    lp=ed_.get(qp,qs);
                    ls=ed_.get(qp,nu[qp]+qs);
                    n_copy(j,i,qp,qs);
                    ed_.set(j,i,lp);
                    ed_.set(j,nu[j]+i,ls);
                    ed_.set(lp,ls,j);
                    ed_.set(lp,nu[lp]+ls,i);
                    ed_.set(qp,qs,-1);
                    i++;
                }
                qs=cycle_up(qs,qp);
                cs=i;
                cp=j;
                n_copy(j,new_double_edge?0:cs,qp,qs);

                // Update the double_edge flag, to pass it
                // to the next instance of this routine
                double_edge=new_double_edge;
            }
        }

        // Connect the final created vertex to the initial one
        ed_.set(cp,cs,rp);
	    ed_.set(rp,0,cp);
        ed_.set(cp,nu[cp]+cs,0);
        ed_.set(rp,nu[rp],cs);
        return false;
    }
    private boolean collapse_order1() {
        while(mec[1]>0) {
            up=0;
            if (Config.VOROPP_VERBOSE >= 1) {
                System.err.println("Order one collapse");
            }
            int i=--mec[1];
            int j=mep[1][3*i];
            int k=mep[1][3*i+1];
            i=mep[1][3*i+2];
            if(!delete_connection(j,k,false)) return false;
            --p;
            if(up==i) up=0;
            if(p!=i) {
                if(up==p) up=i;
                pts[i<<2]=pts[p<<2];
                pts[(i<<2)+1]=pts[(p<<2)+1];
                pts[(i<<2)+2]=pts[(p<<2)+2];
                for(k=0;k<nu[p];k++) ed_.set(ed_.get(p,k), ed_.get(p,nu[p]+k), i);
                n_copy_pointer(i,p);
                ed_.setIndex(i, ed_.storage[p], ed_.getOffset(p));
                nu[i]=nu[p];
                ed_.set(i,nu[i]<<1,i);
            }
        }
        return true;
    }

    /** During the creation of a new facet in the plane routine, it is possible
     * that some order two vertices may arise. This routine removes them.
     * Suppose an order two vertex joins c and d. If there's a edge between
     * c and d already, then the order two vertex is just removed; otherwise,
     * the order two vertex is removed and c and d are joined together directly.
     * It is possible this process will create order two or order one vertices,
     * and the routine is continually run until all of them are removed.
     * \return False if the vertex removal was unsuccessful, indicative of the cell
     *         reducing to zero volume and disappearing; true if the vertex removal
     *         was successful. */
    private boolean collapse_order2() {
        if(!collapse_order1()) return false;
        int a,b,i,j,k,l;
        while(mec[2]>0) {

            // Pick a order 2 vertex and read in its edges
            i=--mec[2];
            j=mep[2][5*i];k=mep[2][5*i+1];
            if(j==k) {
                if (Config.VOROPP_VERBOSE >=1) {
                    System.err.println("Order two vertex joins itself");
                }
                return false;
            }

            // Scan the edges of j to see if joins k
            for(l=0;l<nu[j];l++) {
                if(ed_.get(j,l)==k) break;
            }

            // If j doesn't already join k, join them together.
            // Otherwise delete the connection to the current
            // vertex from j and k.
            a=mep[2][5*i+2];b=mep[2][5*i+3];i=mep[2][5*i+4];
            if(l==nu[j]) {
                ed_.set(j,a,k);
                ed_.set(k,b,j);
                ed_.set(j,nu[j]+a,b);
                ed_.set(k,nu[k]+b,a);
            } else {
                if(!delete_connection(j,a,false)) return false;
                if(!delete_connection(k,b,true)) return false;
            }

            // Compact the memory
            --p;
            if(up==i) up=0;
            if(p!=i) {
                if(up==p) up=i;
                pts[3*i]=pts[3*p];
                pts[3*i+1]=pts[3*p+1];
                pts[3*i+2]=pts[3*p+2];
                for(k=0;k<nu[p];k++) ed_.set(ed_.get(p,k), ed_.get(p,nu[p]+k), i);
                n_copy_pointer(i,p);
                ed_.setIndex(i, ed_.storage[i], ed_.getOffset(p));
                nu[i]=nu[p];
                ed_.set(i,nu[i]<<1,i);
            }

            // Collapse any order 1 vertices if they were created
            if(!collapse_order1()) return false;
        }
        return true;
    }

    /** This routine deletes the kth edge of vertex j and reorganizes the memory.
     * If the neighbor computation is enabled, we also have to supply an handedness
     * flag to decide whether to preserve the plane on the left or right of the
     * connection.
     * \return False if a zero order vertex was formed, indicative of the cell
     *         disappearing; true if the vertex removal was successful. */
    private boolean delete_connection(int j,int k,boolean hand) {
        int q=hand?k:cycle_up(k,j);
        int i=nu[j]-1,l;
        int m;
        if (Config.VOROPP_VERBOSE >= 1) {
            if (i < 1) {
                System.err.println("Zero order vertex formed");
                return false;
            }
        }
        if(mec[i]==mem[i]) add_memory(i);
        n_set_aux1(i);
        for(l=0;l<q;l++) n_copy_aux1(j,l);
        while(l<i) {
            n_copy_aux1_shift(j,l);
            l++;
        }
        Int1Darray edp_ = new Int1Darray(mep[i], ((i<<1)+1)*mec[i]);
        mec[i]++;
        edp_.set(i<<i, j);
        for(l=0;l<k;l++) {
            edp_.set(l, ed_.get(j,l));
            edp_.set(l+i, ed_.get(j,l+nu[j]));
        }
        while(l<i) {
            m=ed_.get(j,l+1);
            edp_.set(l,m);
            k=ed_.get(j,l+nu[j]+1);
            edp_.set(l+i,k);
            ed_.set(m,nu[m]+k, ed_.get(m,nu[m]+k)-1);
            l++;
        }

        mec[nu[j]]--;
        Int1Darray edd_ = new Int1Darray(mep[nu[j]], ((nu[j]<<1)+1)*mec[nu[j]]);
        for(l=0;l<=(nu[j]<<1);l++) ed_.set(j,l,edd_.get(l));
        n_set_aux2_copy(j,nu[j]);
        n_copy_pointer(edd_.get(nu[j]<<1),j);
        n_set_to_aux1(j);
        ed_.setIndex(edd_.get(nu[j]<<1), ed_.getStorage(j), ed_.getOffset(j));
        ed_.setIndex(j, edp_.storage, edp_.getOffset());
        nu[j]=i;
        return true;
    }

    /** Starting from a point within the current cutting plane, this routine attempts
     * to find an edge to a point outside the cutting plane. This prevents the plane
     * routine from .
     * \param[in,out] up */
    private boolean search_for_outside_edge(int[] up) {
        int i,lp,lw;
        int j = stackp2;
        int sc2 = stackp2;
        ds2[stackp2] = up[0];
        stackp2++;
        while(j<stackp2) {
            up[0]=ds2[j];
            j++;
            for(i=0;i<nu[up[0]];i++) {
                lp=ed_.get(up[0], i);
                lw=m_test(lp,new double[1]);
                if (lw==0) {
                    stackp2 = sc2;
                    return true;
                }
                else if (lw==1) add_to_stack(sc2, lp);
            }
        }
        stackp2 = sc2;
        return false;
    }

    /** Adds a point to the auxiliary delete stack if it is not already there.
     * \param[in] vc a reference to the specialized version of the calling class.
     * \param[in] lp the index of the point to add.
     * \param[in,out] stackp2 a pointer to the end of the stack entries. */
    private void add_to_stack(int sc2, int lp) {
        for (int k=sc2; k<stackp2; k++) if (ds2[k] == lp) return;
        if(stackp2==ds2.length) add_memory_ds2();
        ds2[stackp2] = lp;
        stackp2++;
    }
    private void reset_mask() {
        for (int i=0; i<current_vertices; i++) mask[i] = 0;
        maskc = 4;
    }

    /** Assuming that the point up is outside the cutting plane, this routine
     * searches upwards along edges trying to find an edge that intersects the
     * cutting plane.
     * \param[in] rsq the distance along this vector of the plane.
     * \param[in,out] u the dot product of point up with the normal.
     * \return True if the cutting plane was reached, false otherwise. */
    public boolean search_upward(int[] uw,int[] lp,int[] ls,int[] us,double[] l,double[] u) {
        int vs;
        lp[0]=up;l[0]=u[0];

        // The test point is outside of the cutting space
        for(ls[0]=0;ls[0]<nu[lp[0]];ls[0]++) {
               up=ed_.get(lp[0],ls[0]);
               uw[0]=m_test(up,u);
               if(u[0]>l[0]) break;
        }
        if(ls[0]==nu[lp[0]]) if(definite_max(lp,ls,l,u,uw)) {
               up=lp[0];
               return false;
        }

        while(uw[0]==0) {
               //if(++count>=p) failsafe_find(lp,ls,us,l,u);

               // Test all the neighbors of the current point
               // and find the one which is closest to the
               // plane
               vs=ed_.get(lp[0], nu[lp[0]]+ls[0]);
               lp[0]=up;l[0]=u[0];
               for(ls[0]=0;ls[0]<nu[lp[0]];ls[0]++) {
                       if(ls[0]==vs) continue;
                       up=ed_.get(lp[0],ls[0]);
                       uw[0]=m_test(up,u);
                       if(u[0]>l[0]) break;
               }
               if(ls[0]==nu[lp[0]]&&definite_max(lp,ls,l,u,uw)) {
                       up=lp[0];
                       return false;
               }
        }
        us[0]=ed_.get(lp[0], nu[lp[0]]+ls[0]);
        return true;
    }

    /** Checks whether a particular point lp is a definite maximum, searching
     * through any possible minor non-convexities, for a better maximum.
     * \param[in] (x,y,z) the normal vector to the plane. */
    private boolean definite_max(int[] lp,int[] ls,double[] l,double[] u,int[] uw) {
       int tp=lp[0],ts,qp=0;
       int qw;
       double[] q = new double[1];

       // Check to see whether point up is a well-defined maximum. Otherwise
       // any neighboring vertices of up that are marginal need to be
       // followed, to see if they lead to a better maximum.
       for(ts=0;ts<nu[tp];ts++) {
               qp=ed_.get(tp,ts);
               m_test(qp,q);
               if(q[0]>l[0]-big_tol) break;
       }
       if(ts==nu[tp]) return true;

       // The point tp is marginal, so it will be necessary to do the
       // flood-fill search. Mark the point tp and the point qp, and search
       // any remaining neighbors of the point tp.
       int stackp=1;
       flip(lp[0]);
       flip(qp);
       ds[0]=qp;
       ts++;
       while(ts<nu[tp]) {
               qp=ed_.get(tp,ts);
               m_test(qp,q);
               if(q[0]>l[0]-big_tol) {
                       if(stackp==ds.length) add_memory_ds();
                       ds[stackp] = up;
                       stackp++;
                       flip(up);
               }
               ts++;
       }

       // Consider additional marginal points, starting with the original
       // point qp
       int spp=0;
       while(spp<stackp) {
           tp = ds[spp];
           spp++;
               for(ts=0;ts<nu[tp];ts++) {
                       qp=ed_.get(tp,ts);

                       // Skip the point if it's already marked
                       if(ed_.get(qp,nu[qp]<<1)<0) continue;
                       qw=m_test(qp,q);

                       // This point is a better maximum. Reset markers and
                       // return true.
                       if(q[0]>l[0]) {
                               flip(lp[0]);
                               lp[0]=tp;
                               ls[0]=ts;
                               m_test(lp[0],l);
                               up=qp;
                               uw[0]=qw;
                               u=q;
                               while(stackp>0) {
                                   stackp--;
                                   flip(ds[stackp]);
                               };
                               return false;
                       }

                       // The point is marginal and therefore must also be
                       // considered
                       if(q[0]>l[0]-big_tol) {
                               if(stackp==ds.length) {
                                       int nn=stackp-spp;
                                       add_memory_ds();
                                       spp=stackp-nn;
                               }
                               ds[stackp] = qp;
                               stackp++;
                               flip(qp);
                       }
               }
       }

       // Reset markers and return false
       flip(lp[0]);
       while(stackp>0) {
           stackp--;
           flip(ds[stackp]);
       }
       return true;
    }

    private boolean search_downward(int[] lw,int[] lp,int[] ls,int[] us,double[] l,double[] u) {
       int vs;

       // The test point is outside of the cutting space
       for(us[0]=0;us[0]<nu[up];us[0]++) {
               lp[0]=ed_.get(up,us[0]);
               lw[0]=m_test(lp[0],l);
               if(u[0]>l[0]) break;
       }
       if(us[0]==nu[up]) if(definite_min(lp,us,l,u,lw)) return false;

       while(lw[0]==2) {
               //if(++count>=p) failsafe_find(lp,ls,us,l,u);

               // Test all the neighbors of the current point
               // and find the one which is closest to the
               // plane
               vs=ed_.get(up,nu[up]+us[0]);
               up=lp[0];
               u[0]=l[0];
               for(us[0]=0;us[0]<nu[up];us[0]++) {
                       if(us[0]==vs) continue;
                       lp[0]=ed_.get(up,us[0]);
                       lw[0]=m_test(lp[0],l);
                       if(u[0]>l[0]) break;
               }
               if(us[0]==nu[up]&&definite_min(lp,us,l,u,lw)) return false;
       }
       ls[0]=ed_.get(up,nu[up]+us[0]);
       return true;
    }
    private boolean definite_min(int[] lp,int[] us,double[] l,double[] u,int[] lw) {
        int tp=up,ts,qp=0;
        int qw;
        double[] q = new double[1];

        // Check to see whether point up is a well-defined maximum. Otherwise
        // any neighboring vertices of up that are marginal need to be
       // followed, to see if they lead to a better maximum.
        for(ts=0;ts<nu[tp];ts++) {
                qp=ed_.get(tp,ts);
                m_test(qp,q);
                if(q[0]<u[0]+big_tol) break;
        }
        if(ts==nu[tp]) return true;

        // The point tp is marginal, so it will be necessary to do the
        // flood-fill search. Mark the point tp and the point qp, and search
        // any remaining neighbors of the point tp.
        int stackp=1;
        flip(up);
        flip(qp);
        ds[0]=qp;
        ts++;
        while(ts<nu[tp]) {
                qp=ed_.get(tp,ts);
                m_test(qp,q);
                if(q[0]<u[0]+big_tol) {
                        if(stackp==ds.length) add_memory_ds();
                        ds[stackp] = lp[0];
                        stackp++;
                        flip(lp[0]);
                }
                ts++;
        }

        // Consider additional marginal points, starting with the original
        // point qp
        int spp=0; //ds;
        while(spp<stackp) {
                tp=ds[spp];
                spp++;
                for(ts=0;ts<nu[tp];ts++) {
                        qp=ed_.get(tp,ts);

                        // Skip the point if it's already marked
                        if(ed_.get(qp,nu[qp]<<1)<0) continue;
                        qw=m_test(qp,q);

                        // This point is a better minimum. Reset markers and
                        // return true.
                        if(q[0]<u[0]) {
                                flip(up);
                                up=tp;
                                us[0]=ts;
                                m_test(up,u);
                                lp[0]=qp;
                                lw[0]=qw;
                                l=q;
                                while(stackp>0) {
                                    stackp--;
                                    flip(ds[stackp]);
                                };
                                return false;
                        }

                        // The point is marginal and therefore must also be
                        // considered
                        if(q[0]<u[0]+big_tol) {
                                if(stackp==ds.length) {
                                        int nn=stackp-spp;
                                        add_memory_ds();
                                        spp=stackp-nn;
                                }
                                ds[stackp] = qp;
                                stackp++;
                                flip(qp);
                        }
                }
        }

        // Reset markers and return false
        flip(up);
        while(stackp>0) {
            stackp--;
            flip(ds[stackp]);
        };
        return true;
    }
    private void minkowski_contrib(int i,int k,int m,double r,double[] ar,double[] vo) {
        double ix=pts[4*i],iy=pts[4*i+1],iz=pts[4*i+2],
                kx=pts[4*k],ky=pts[4*k+1],kz=pts[4*k+2],
                mx=pts[4*m],my=pts[4*m+1],mz=pts[4*m+2],
                ux=kx-ix,uy=ky-iy,uz=kz-iz,vx=mx-kx,vy=my-ky,vz=mz-kz,
                e1x=uz*vy-uy*vz,e1y=ux*vz-uz*vx,e1z=uy*vx-ux*vy,e2x,e2y,e2z,
                wmag=e1x*e1x+e1y*e1y+e1z*e1z;
        if(wmag<tol*tol) return;
        wmag=1/Math.sqrt(wmag);
        e1x*=wmag;e1y*=wmag;e1z*=wmag;

        // Compute second orthonormal vector
        if(Math.abs(e1x)>0.5) {
            e2x=-e1y;e2y=e1x;e2z=0;
        } else if(Math.abs(e1y)>0.5) {
            e2x=0;e2y=-e1z;e2z=e1y;
        } else {
            e2x=e1z;e2y=0;e2z=-e1x;
        }
        wmag=1/Math.sqrt(e2x*e2x+e2y*e2y+e2z*e2z);
        e2x*=wmag;e2y*=wmag;e2z*=wmag;

        // Compute third orthonormal vector
        double e3x=e1z*e2y-e1y*e2z,
                e3y=e1x*e2z-e1z*e2x,
                e3z=e1y*e2x-e1x*e2y,
                x0=e1x*ix+e1y*iy+e1z*iz;
        if(x0<tol) return;

        double ir=e2x*ix+e2y*iy+e2z*iz,is=e3x*ix+e3y*iy+e3z*iz,
                kr=e2x*kx+e2y*ky+e2z*kz,ks=e3x*kx+e3y*ky+e3z*kz,
                mr=e2x*mx+e2y*my+e2z*mz,ms=e3x*mx+e3y*my+e3z*mz;

        minkowski_edge(x0,ir,is,kr,ks,r,ar,vo);
        minkowski_edge(x0,kr,ks,mr,ms,r,ar,vo);
        minkowski_edge(x0,mr,ms,ir,is,r,ar,vo);
    }
    private void minkowski_edge(double x0,double r1,double s1,double r2,double s2,double r,double[] ar,double[] vo) {
        double r12=r2-r1,s12=s2-s1,l12=r12*r12+s12*s12;
        if(l12<tol*tol) return;
        l12=1/Math.sqrt(l12);r12*=l12;s12*=l12;
        double y0=s12*r1-r12*s1;
        if(Math.abs(y0)<tol) return;
        minkowski_formula(x0,y0,-r12*r1-s12*s1,r,ar,vo);
        minkowski_formula(x0,y0,r12*r2+s12*s2,r,ar,vo);
    }
    private void minkowski_formula(double x0,double y0,double z0,double r,double[] ar,double[] vo) {
	    final double pi=3.1415926535897932384626433832795;
        if(Math.abs(z0)<tol) return;
        double si;
        if(z0<0) {z0=-z0;si=-1;} else si=1;
        if(y0<0) {y0=-y0;si=-si;}
        double xs=x0*x0,ys=y0*y0,zs=z0*z0,res=xs+ys,rvs=res+zs,theta=Math.atan(z0/y0),rs=r*r,rc=rs*r,temp,voc,arc;
        if(r<x0) {
            temp=2*theta-0.5*pi-Math.asin((zs*xs-ys*rvs)/(res*(ys+zs)));
            voc=rc/6.*temp;
            arc=rs*0.5*temp;
        } else if(rs<res*1.0000000001) {
            temp=0.5*pi+Math.asin((zs*xs-ys*rvs)/(res*(ys+zs)));
            voc=theta*0.5*(rs*x0-xs*x0/3.)-rc/6.*temp;
            arc=theta*x0*r-rs*0.5*temp;
        } else if(rs<rvs) {
            temp=theta-pi*0.5+Math.asin(y0/Math.sqrt(rs-xs));
            double temp2=(rs*x0-xs*x0/3.),
                    x2s=rs*xs/res,y2s=rs*ys/res,
                    temp3=Math.asin((x2s-y2s-xs)/(rs-xs)),
                    temp4=Math.asin((zs*xs-ys*rvs)/(res*(ys+zs))),
                    temp5=Math.sqrt(rs-res);
            voc=0.5*temp*temp2+x0*y0/6.*temp5+r*rs/6*(temp3-temp4);
            arc=x0*r*temp-0.5*temp2*y0*r/((rs-xs)*temp5)+x0*y0/6.*r/temp5+rs*0.5*temp3+rs*rs/3.*2*xs*ys/(res*(rs-xs)*Math.sqrt((rs-xs)*(rs-xs)-(x2s-y2s-xs)*(x2s-y2s-xs)))-rs*0.5*temp4;
        } else {
            voc=x0*y0*z0/6.;
            arc=0;
        }
        vo[0]+=voc*si;
        ar[0]+=arc*si;
    }

    public static double dot_product(double[] a, int aoff, double[] b, int boff) {
        return a[aoff+0]*b[boff+0]+a[aoff+1]*b[boff+1]+a[aoff+2]*b[boff+2];
    }

    public static void cross_product(double[] a, double[] b, double[] c) {
        c[0]=a[1]*b[2]-a[2]*b[1];
        c[1]=a[2]*b[0]-a[0]*b[2];
        c[2]=a[0]*b[1]-a[1]*b[0];
    }

    public static double triple_product(double[] a, double[] b, double[] c) {
        double[] cp = new double[3];
        cross_product(b,c,cp);
        return dot_product(a,0,cp,0);
    }

    public static void normalize_vector(double[] x, int offset, double[] normalized) {
        double normsq = dot_product(x,offset,x,offset);
        double invnorm = normsq<=0?0:1/Math.sqrt(normsq);
        for (int i=0;i<3;i++) {
            normalized[i]=x[i]*invnorm;
        }
    }

    public static double calculate_solid_angle(double[] R1, double[] R2, double[] R3) {
        // Method described in:
        // A. Van Oosterom and J. Strackee
        // "The Solid Angle of a Plane Triangle"
        // IEEE Transactions on Biomedical Engineering, BME-30, 2, 1983, 125--126
        // https://doi.org/10.1109/TBME.1983.325207
        return Math.abs(2*Math.atan2(triple_product(R1,R2,R3),1+dot_product(R1,0,R2,0)+dot_product(R2,0,R3,0)+dot_product(R3,0,R1,0)));
    }

    /* This routine tests to see if a cell intersects a plane, by tracing over the
     * cell from vertex to vertex, starting at up. It is meant to be called either
     * by plane_intersects() or plane_intersects_track(), when those routines
     * cannot immediately resolve the case.
     * \param[in] (x,y,z) the normal vector to the plane.
     * \param[in] rsq the distance along this vector of the plane.
     * \param[in] g the distance of up from the plane.
     * \return False if the plane does not intersect the plane, true if it does. */
    private boolean plane_intersects_track(double x,double y,double z,double rsq,double g) {
        for(int tp=0;tp<p;tp++) if(x*pts[tp<<2]+y*pts[(tp<<2)+1]+z*pts[(tp<<2)+2]>rsq) return true;
        return false;
/*
	int ls,us,lp;
	double l,u;
	unsigned int uw;

	// Initialize the safe testing routine
	px=x;py=y;pz=z;prsq=rsq;
	maskc+=4;
	if(maskc<4) reset_mask();

	return search_upward(uw,lp,ls,us,l,u);
}*/

        /*
        int count=0,ls,us,tp;
        double t;

        // The test point is outside of the cutting space
        for(us=0;us<nu[up];us++) {
            tp=ed_.get(up,us);
            t=x*pts[4*tp]+y*pts[4*tp+1]+z*pts[4*tp+2];
            if(t>g) {
                ls=ed_.get(up,nu[up]+us);
                up=tp;
                while (t<rsq) {
                    if(++count>=p) {
                        if (Config.VOROPP_VERBOSE >= 1) {
                            System.err.println("Bailed out of convex calculation");
                        }
                        for(tp=0;tp<p;tp++) if(x*pts[4*tp]+y*pts[4*tp+1]+z*pts[4*tp+2]>rsq) return true;
                        return false;
                    }

                    // Test all the neighbors of the current point
                    // and find the one which is closest to the
                    // plane
                    for(us=0;us<ls;us++) {
                        tp=ed_.get(up,us);
                        g=x*pts[4*tp]+y*pts[4*tp+1]+z*pts[4*tp+2];
                        if(g>t) break;
                    }
                    if(us==ls) {
                        us++;
                        while(us<nu[up]) {
                            tp=ed_.get(up,us);
                            g=x*pts[4*tp]+y*pts[4*tp+1]+z*pts[4*tp+2];
                            if(g>t) break;
                            us++;
                        }
                        if(us==nu[up]) return false;
                    }
                    ls=ed_.get(up,nu[up]+us);up=tp;t=g;
                }
                return true;
            }
        }
        return false;*/
    }

    /** This inline routine is called by normals(). It attempts to construct a
     * single normal vector that is associated with a particular face. It first
     * traces around the face, trying to find two vectors along the face edges
     * whose vector product is above the numerical tolerance. It then constructs
     * the normal vector using this product. If the face is too small, and none of
     * the vector products are large enough, the routine may return (0,0,0) as the
     * normal vector.
     * \param[in] v the vector to store the results in.
     * \param[in] i the initial vertex of the face to test.
     * \param[in] j the index of an edge of the vertex.
     * \param[in] k the neighboring vertex of i, set to ed[i][j]. */
    private void normals_search(DoubleArrayList v,int i,int j,int k) {
        ed_.set(i,j,-1-k);;
        int l=cycle_up(ed_.get(i,nu[i]+j),k),m;
        double ux,uy,uz,vx,vy,vz,wx,wy,wz,wmag;
        do {
            m=ed_.get(k,l);ed_.set(k,l,-1-m);
            ux=pts[4*m]-pts[4*k];
            uy=pts[4*m+1]-pts[4*k+1];
            uz=pts[4*m+2]-pts[4*k+2];

            // Test to see if the length of this edge is above the tolerance
            if(ux*ux+uy*uy+uz*uz>tol) {
                while(m!=i) {
                    l=cycle_up(ed_.get(k,nu[k]+l),m);
                    k=m;m=ed_.get(k,l);ed_.set(k,l,-1-m);
                    vx=pts[4*m]-pts[4*k];
                    vy=pts[4*m+1]-pts[4*k+1];
                    vz=pts[4*m+2]-pts[4*k+2];

                    // Construct the vector product of this edge with
                    // the previous one
                    wx=uz*vy-uy*vz;
                    wy=ux*vz-uz*vx;
                    wz=uy*vx-ux*vy;
                     wmag=wx*wx+wy*wy+wz*wz;

                    // Test to see if this vector product of the
                    // two edges is above the tolerance
                    if(wmag>tol) {

                        // Construct the normal vector and print it
                        wmag=1/Math.sqrt(wmag);
                        v.add(wx*wmag);
                        v.add(wy*wmag);
                        v.add(wz*wmag);

                        // Mark all of the remaining edges of this
                        // face and exit
                        while(m!=i) {
                            l=cycle_up(ed_.get(k,nu[k]+l),m);
                            k=m;m=ed_.get(k,l);ed_.set(k,l,-1-m);
                        }
                        return;
                    }
                }
                v.add(0);
                v.add(0);
                v.add(0);
                return;
            }
            l=cycle_up(ed_.get(k,nu[k]+l),m);
            k=m;
        } while (k!=i);
        v.add(0);
        v.add(0);
        v.add(0);
    }
    private int search_edge(int l) {
        for(int m=0;m<nu[l];m++) {
            int k=ed_.get(l,m);
            if(k>=0) return m;
        }
        return -1;
    }

    /** Checks to see if a given vertex is inside, outside or within the test
     * plane. If the point is far away from the test plane, the routine immediately
     * returns whether it is inside or outside. If the routine is close the the
     * plane and within the specified tolerance, then the special check_marginal()
     * routine is called.
     * \param[in] n the vertex to test.
     * \param[out] ans the result of the scalar product used in evaluating the
     *                 location of the point.
     * \return -1 if the point is inside the plane, 1 if the point is outside the
     *         plane, or 0 if the point is within the plane. */
    private int m_test(int n,double[] ans) {
        if(mask[n]>=maskc) {
                ans[0]=pts[4*n+3];
                return mask[n]&3;
        } else return m_calc(n,ans);
    }

    /** Checks to see if a given vertex is inside, outside or within the test
     * plane. If the point is far away from the test plane, the routine immediately
     * returns whether it is inside or outside. If the routine is close the the
     * plane and within the specified tolerance, then the special check_marginal()
     * routine is called.
     * \param[in] n the vertex to test.
     * \param[out] ans the result of the scalar product used in evaluating the
     *                 location of the point.
     * \return -1 if the point is inside the plane, 1 if the point is outside the
     *         plane, or 0 if the point is within the plane. */
    private int m_testx(int n, double[] ans) {
        int maskr;
        if(mask[n]>=maskc) {
                ans[0]=pts[4*n+3];
                maskr=mask[n]&3;
        } else maskr=m_calc(n,ans);
        if(maskr==0&&ans[0]>-big_tol&&ed_.get(n,nu[n]<<1)!=-1) {
            ed_.set(n, nu[n]<<1, -1);
            if(stackp3==xse.length) add_memory_xse();
            xse[stackp3] = n;
            stackp3++;
        }
        return maskr;
    }
    private int m_calc(int n, double[] ans) {
        int idx = 4*n;
        ans[0]=pts[idx+0]*px + pts[idx+1]*py + pts[idx+2]*pz - prsq;
        pts[idx+3] = ans[0];
        int maskr = ans[0] < -tol ? 0 : (ans[0]>tol ? 2 : 1);
        mask[n] = maskc | maskr;
        return maskr;
    }
    private void flip(int tp) {
        ed_.set(tp, nu[tp]<<1, -1-ed_.get(tp,nu[tp]<<1));
    }
}
