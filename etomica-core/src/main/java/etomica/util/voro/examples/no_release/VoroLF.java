// File import example code
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro.examples.no_release;

import etomica.util.collections.DoubleArrayList;
import etomica.util.collections.IntArrayList;
import etomica.util.voro.*;
import etomica.util.voro.examples.ResourceHelper;

import java.io.*;

public class VoroLF {
    
    public static void main(String[] args) {
        // Set up constants for the container geometry
        final double ax=-0.5,bx=25.5;
        final double ay=-0.5,by=25.5;
        final double az=-0.5,bz=25.5;

        // Manually import the file
        int i,j,max_id=0,n;
        IntArrayList vid = new IntArrayList();
        IntArrayList neigh = new IntArrayList();
        IntArrayList f_order = new IntArrayList();
        DoubleArrayList vx = new DoubleArrayList();
        DoubleArrayList vy = new DoubleArrayList();
        DoubleArrayList vz = new DoubleArrayList();
        DoubleArrayList vd = new DoubleArrayList();
        try {
            InputStream fp = ResourceHelper.getStreamForFile("liq-900K.dat", VoroLF.class);
            BufferedReader bufReader = new BufferedReader(new InputStreamReader(fp));

            String line = null;
            while ((line = bufReader.readLine()) != null) {
                String[] bits = line.trim().split("[ \t]+");
                if (bits.length != 4) Common.voro_fatal_error("File import error", Config.Voropp.FILE_ERROR);
                int id = Integer.parseInt(bits[0]);
                vid.add(id);
                max_id = Math.max(max_id, id);
                double x = Double.parseDouble(bits[1]);
                double y = Double.parseDouble(bits[2]);
                double z = Double.parseDouble(bits[3]);
                vx.add(x);
                vy.add(y);
                vz.add(z);
            }
            fp.close();
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
        n = vid.size();

        // Compute optimal size for container, and then construct the container
        double dx=bx-ax,dy=by-ay,dz=bz-az;
        double l = Math.pow(n/(5.6*dx*dy*dz),1/3.0);
        int nx=(int)(dx*l+1);
        int ny=(int)(dy*l+1);
        int nz=(int)(dz*l+1);
        Container con = new Container(ax,bx,ay,by,az,bz,nx,ny,nz,false,false,false,8);

        // Print status message
        System.out.printf("Read %d particles, max ID is %d\n"+
                "Container grid is %d by %d by %d\n",n,max_id,nx,ny,nz);

        // Import the particles, and create ID lookup tables
        double[] xi=new double[max_id+1];
        double[] yi=new double[max_id+1];
        double[] zi=new double[max_id+1];
        for(j=0;j<n;j++) {
            int id=vid.getInt(j);
            double x=vx.getDouble(j);
            double y=vy.getDouble(j);
            double z=vz.getDouble(j);
            con.put(id,x,y,z);
            xi[id]=x;
            yi[id]=y;
            zi[id]=z;
        }

        try {
            // Open three output files for statistics and gnuplot cells
            FileOutputStream fp = new FileOutputStream("liq-900K.out");
            FileOutputStream fp2 = new FileOutputStream("liq-900K.gnu");
            FileOutputStream fp3 = new FileOutputStream("liq-900K-orig.gnu");

            // Loop over all particles and compute their Voronoi cells
            VoronoiCellNeighbor c = new VoronoiCellNeighbor();
            VoronoiCellNeighbor c2 = new VoronoiCellNeighbor();
            CLoopAll cl = new CLoopAll(con);
            if (cl.start()) do if (con.compute_cell(c, cl)) {

                // Get particle position, ID, and neighbor vector
                double[] xyz = cl.pos();
                double x = xyz[0];
                double y = xyz[1];
                double z = xyz[2];
                int id = cl.pid();
                c.neighbors(neigh);

                // Get face areas et total surface of faces
                c.face_areas(vd);
                c.surface_area();
                c.draw_gnuplot(x, y, z, fp3);

                // Initialize second cell
                c2.init(ax - x, bx - x, ay - y, by - y, az - z, bz - z);

                // Add condition on surface: >1% total surface. In addition,
                // skip negative indices, since they correspond to faces
                // against the container boundaries
                for (i = 0; i < vd.size(); i++) {
                    if (vd.getDouble(i) > 0.01 * c.surface_area() && neigh.getInt(i) >= 0) {
                        j = neigh.getInt(i);
                        c2.nplane(xi[j] - x, yi[j] - y, zi[j] - z, j);
                    }
                }

                // Get information of c2 cell
                c2.face_areas(vd);
                c2.face_orders(f_order);

                // Output information to file
                i = vd.size();
                fp.write(String.format("%d %d", id, i).getBytes());
                for (j = 0; j < i; j++) fp.write(String.format(" %d", f_order.getInt(j)).getBytes());
                for (j = 0; j < i; j++) fp.write(String.format(" %.3f", vd.getDouble(j)).getBytes());
                fp.write(String.format(" %.3f %.3f %.3f\n", x, y, z).getBytes());

                c2.draw_gnuplot(x, y, z, fp2);
            } while (cl.inc());

            // Close files
            fp.close();
            fp2.close();
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }

    }
}
