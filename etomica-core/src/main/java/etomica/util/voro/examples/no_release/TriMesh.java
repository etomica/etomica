// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro.examples.no_release;

import etomica.util.collections.DoubleArrayList;
import etomica.util.collections.IntArrayList;
import etomica.util.voro.CLoopAll;
import etomica.util.voro.Container;
import etomica.util.voro.VoronoiCell;

public class TriMesh {

    public static void main(String[] args) {

        // Set up constants for the container geometry
        final double x_min=-1,x_max=1;
        final double y_min=-1,y_max=1;
        final double z_min=-1,z_max=1;
        
        // Set up the number of blocks that the container is divided into
        final int n_x=6,n_y=6,n_z=6;
        
        // Set the number of particles that are going to be randomly introduced
        final int particles=1;

        int i,j,id;

        // Create a container with the geometry given above, and make it
        // non-periodic in each of the three coordinates. Allocate space for
        // eight particles within each computational block
        Container con = new Container(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
                false,false,false,8);

        // Randomly add particles into the container
        for(i=0;i<particles;i++) {
            double x=x_min+Math.random()*(x_max-x_min);
            double y=y_min+Math.random()*(y_max-y_min);
            double z=z_min+Math.random()*(z_max-z_min);
            con.put(i,x,y,z);
        }

        // Sum up the volumes, and check that this matches the container volume
        CLoopAll cl = new CLoopAll(con);
        IntArrayList f_vert = new IntArrayList();
        DoubleArrayList v = new DoubleArrayList();
        VoronoiCell c = new VoronoiCell();
        if(cl.start()) do if(con.compute_cell(c,cl)) {
            double[] x = new double[1];
            double[] y = new double[1];
            double[] z = new double[1];
            cl.pos(x,y,z);
            id=cl.pid();
            System.out.printf("Particle %d:\n",id);

            // Gather information about the computed Voronoi cell
            c.face_vertices(f_vert);
            c.vertices(x[0],y[0],z[0],v);

            // Print vertex positions
            for(i=0;i<v.size();i+=3) System.out.printf("Vertex %d : (%g,%g,%g)\n",i/3,v.getDouble(i),v.getDouble(i+1),v.getDouble(i+2));
            System.out.println();

            // Loop over all faces of the Voronoi cell
            j=0;
            while(j<f_vert.size()) {

                // Number of vertices in this face
                int nv=f_vert.getInt(j);

                // Print triangles
                for(i=2;i<nv;i++)
                    System.out.printf("Triangle : (%d,%d,%d)\n",f_vert.getInt(j+1),f_vert.getInt(j+i),f_vert.getInt(j+i+1));

                // Move j to point at the next face
                j+=nv+1;
            }
            System.out.println();
        } while (cl.inc());

        // Output the particle positions in gnuplot format
        con.draw_particles("random_points_p.gnu");

        // Output the Voronoi cells in gnuplot format
        con.draw_cells_gnuplot("random_points_v.gnu");

    }
}
