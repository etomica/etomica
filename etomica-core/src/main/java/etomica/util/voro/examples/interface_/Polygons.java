// Direct C++ interface example code
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro.examples.interface_;

import etomica.util.collections.DoubleArrayList;
import etomica.util.collections.IntArrayList;
import etomica.util.voro.CLoopAll;
import etomica.util.voro.Container;
import etomica.util.voro.PreContainer;
import etomica.util.voro.VoronoiCellNeighbor;
import etomica.util.voro.examples.ResourceHelper;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;

public class Polygons {

    public static void main(String[] args) {

        int i,j;
        int id;
        int[] nx = new int[1];
        int[] ny = new int[1];
        int[] nz = new int[1];
        VoronoiCellNeighbor c = new VoronoiCellNeighbor();
        IntArrayList neigh = new IntArrayList();
        IntArrayList f_vert = new IntArrayList();
        DoubleArrayList v = new DoubleArrayList();

        // Create a pre-container class to import the input file and guess the
        // best computational grid size to use.
        PreContainer pcon = new PreContainer(-3,3,-3,3,0,6,false,false,false);
        InputStream in = ResourceHelper.getStreamForFile("pack_six_cube", Polygons.class);
        pcon.import_(in);
        pcon.guess_optimal(nx,ny,nz);

        // Set up the container class and import the particles from the
        // pre-container
        Container con = new Container(-3,3,-3,3,0,6,nx[0],ny[0],nz[0],false,false,false,8);
        pcon.setup(con);

        // Open the output files
        try {
            FileOutputStream fp4 = new FileOutputStream("polygons4_v.pov");
            FileOutputStream fp5 = new FileOutputStream("polygons5_v.pov");
            FileOutputStream fp6 = new FileOutputStream("polygons6_v.pov");

            // Loop over all particles in the container and compute each Voronoi
            // cell
            CLoopAll cl = new CLoopAll(con);
            if (cl.start()) do if (con.compute_cell(c, cl)) {
                double[] x = new double[1];
                double[] y = new double[1];
                double[] z = new double[1];
                cl.pos(x, y, z);
                id = cl.pid();

                // Gather information about the computed Voronoi cell
                c.neighbors(neigh);
                c.face_vertices(f_vert);
                c.vertices(x[0], y[0], z[0], v);

                // Loop over all faces of the Voronoi cell
                for (i = 0, j = 0; i < neigh.size(); i++) {

                    // Draw all quadrilaterals, pentagons, and hexagons.
                    // Skip if the neighbor information is smaller than
                    // this particle's ID, to avoid double counting. This
                    // also removes faces that touch the walls, since the
                    // neighbor information is set to negative numbers for
                    // these cases.
                    if (neigh.getInt(i) > id) {
                        switch (f_vert.getInt(j)) {
                            case 4:
                                draw_polygon(fp4, f_vert, v, j);
                                break;
                            case 5:
                                draw_polygon(fp5, f_vert, v, j);
                                break;
                            case 6:
                                draw_polygon(fp6, f_vert, v, j);
                        }
                    }

                    // Skip to the next entry in the face vertex list
                    j += f_vert.getInt(j) + 1;
                }
            } while (cl.inc());

            // Close the output files
            fp4.close();
            fp5.close();
            fp6.close();

        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }

        // Draw the outline of the domain
        con.draw_domain_pov("polygons_d.pov");

    }

    public static void draw_polygon(FileOutputStream fp,IntArrayList f_vert,DoubleArrayList v,int j) throws IOException {
        String[] s = new String[6];
        int k,l,n=f_vert.getInt(j);

        // Create POV-Ray vector strings for each of the vertices
        for(k=0;k<n;k++) {
            l=3*f_vert.getInt(j+k+1);
            s[k] = String.format("<%g,%g,%g>",v.getDouble(l),v.getDouble(l+1),v.getDouble(l+2));
        }

        // Draw the interior of the polygon
        fp.write("union{\n".getBytes());
        for(k=2;k<n;k++) fp.write(String.format("\ttriangle{%s,%s,%s}\n",s[0],s[k-1],s[k]).getBytes());
        fp.write("\ttexture{t1}\n}\n".getBytes());

        // Draw the outline of the polygon
        fp.write("union{\n".getBytes());
        for(k=0;k<n;k++) {
            l=(k+1)%n;
            fp.write(String.format("\tcylinder{%s,%s,r}\n\tsphere{%s,r}\n",
                    s[k],s[l],s[l]).getBytes());
        }
        fp.write("\ttexture{t2}\n}\n".getBytes());
    }

}
