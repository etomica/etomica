// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro.examples.no_release;

import etomica.util.collections.IntArrayList;
import etomica.util.voro.CLoopAll;
import etomica.util.voro.Container;
import etomica.util.voro.VoronoiCell;
import etomica.util.voro.VoronoiCellNeighbor;

import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;

public class LlyodBox {
    
    public static void main(String[] args) {
        // Set up constants for the container geometry
        final double boxl=1;

        // Set up the number of blocks that the container is divided into
        final int bl=10;

        // Set the number of particles that are going to be randomly introduced
        final int particles=4000;

        // Set the number of Voronoi faces to bin
        final int nface=40;

        int l;
        int[] faces = new int[nface];
        double[] p = new double[3*particles];

        // Create a container with the geometry given above, and make it
        // non-periodic in each of the three coordinates. Allocate space for
        // eight particles within each computational block
        Container con = new Container(-boxl,boxl,-boxl,boxl,-boxl,boxl,bl,bl,bl,false,false,false,8);

        // Randomly add particles into the container
        for(int i=0;i<particles;i++) {
            double x=boxl*(2*Math.random()-1);
            double y=boxl*(2*Math.random()-1);
            double z=boxl*(2*Math.random()-1);
            con.put(i,x,y,z);
        }

        for(l=0;l<=200;l++) {
            CLoopAll vl = new CLoopAll(con);
            VoronoiCell c = new VoronoiCell();
            Arrays.fill(faces, 0);
            if(vl.start()) do if(con.compute_cell(c,vl)) {
                int[] pid = new int[1];
                double[] x = new double[1];
                double[] y = new double[1];
                double[] z = new double[1];
                double[] r = new double[1];
                vl.pos(pid,x,y,z,r);
                int i = pid[0];
                double[] dxyz = c.centroid();
                p[3*i]=x[0]+dxyz[0];
                p[3*i+1]=y[0]+dxyz[1];
                p[3*i+2]=z[0]+dxyz[2];

                i=c.number_of_faces()-4;
                if(i<0) i=0;if(i>=nface) i=nface-1;
                faces[i]++;
            } while (vl.inc());
            con.clear();
            for(int i=0;i<particles;i++) con.put(i,p[3*i],p[3*i+1],p[3*i+2]);
            System.out.printf("%d",l);
            for(int fp=0;fp<nface;fp++) System.out.printf(" %d",faces[fp]);
            System.out.println();
        }

        // Output the particle positions in gnuplot format
        con.draw_particles("sphere_mesh_p.gnu");

        // Output the Voronoi cells in gnuplot format
        con.draw_cells_gnuplot("sphere_mesh_v.gnu");

        try {
            // Output the neighbor mesh in gnuplot format
            FileOutputStream ff = new FileOutputStream("sphere_mesh.net");
            IntArrayList vi = new IntArrayList();
            VoronoiCellNeighbor c = new VoronoiCellNeighbor();
            CLoopAll vl = new CLoopAll(con);
            if(vl.start()) do if(con.compute_cell(c,vl)) {
                int i=vl.pid();
                c.neighbors(vi);
                for(l=0;l<vi.size();l++) if(vi.getInt(l)>i)
                    ff.write(String.format("%g %g %g\n%g %g %g\n\n\n",
                            p[3*i],p[3*i+1],p[3*i+2],
                            p[3*vi.getInt(l)],p[3*vi.getInt(l)+1],p[3*vi.getInt(l)+2]).getBytes());
            } while (vl.inc());
            ff.close();
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }


    }
}
