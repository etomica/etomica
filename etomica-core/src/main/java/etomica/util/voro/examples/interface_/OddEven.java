// Odd/even face coloring code
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro.examples.interface_;

import etomica.util.collections.DoubleArrayList;
import etomica.util.collections.IntArrayList;
import etomica.util.voro.VoronoiCell;

import java.io.FileOutputStream;
import java.io.IOException;

public class OddEven {
    
    public static void main(String[] args) {
        int i;
        double x,y,z,rsq,r;
        VoronoiCell v = new VoronoiCell();

        // Initialize the Voronoi cell to be a cube of side length 2, centered
        // on the origin
        v.init(-1,1,-1,1,-1,1);

        // Cut the cell by 250 random planes which are all a distance 1 away
        // from the origin, to make an approximation to a sphere
        for(i=0;i<250;i++) {
            x=2*Math.random()-1;
            y=2*Math.random()-1;
            z=2*Math.random()-1;
            rsq=x*x+y*y+z*z;
            if(rsq>0.01&&rsq<1) {
                r=1/Math.sqrt(rsq);x*=r;y*=r;z*=r;
                v.plane(x,y,z,1);
            }
        }

        // Calculate the orders of the faces and the normal vectors
        IntArrayList f_vert = new IntArrayList();
        DoubleArrayList nor = new DoubleArrayList();
        v.face_orders(f_vert);
        v.normals(nor);

        // Output POV-Ray planes with textures based on whether a face is
        // composed of an odd or even number of edges
	    final String[] parity = new String[]{"even","odd"};
        try {
            FileOutputStream fp = new FileOutputStream("odd_even_pl.pov");
            for (i = 0; i < f_vert.size(); i++)
                fp.write(String.format("plane{<%g,%g,%g>,0.5 texture{t_%s}}\n"
                        , nor.getDouble(3 * i), nor.getDouble(3 * i + 1), nor.getDouble(3 * i + 2)
                        , parity[f_vert.getInt(i) & 1]).getBytes());
            fp.close();
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }

        // Save the Voronoi cell as a spheres and cylinders
        v.draw_pov(0,0,0,"odd_even_v.pov");

    }
}
