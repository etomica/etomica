// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro;

import etomica.util.collections.DoubleArrayList;
import etomica.util.collections.IntArrayList;

import java.io.IOException;
import java.io.OutputStream;

public class Common {

    public static void check_duplicate(int n,double x,double y,double z,int id,double[] q, int qp) {
        double dx=q[qp]-x,dy=q[qp+1]-y,dz=q[qp+2]-z;
        if(dx*dx+dy*dy+dz*dz<1e-10) {
            throw new RuntimeException(String.format("Duplicate: %d (%g,%g,%g) matches %d (%g,%g,%g)\n",n,x,y,z,id,q[qp],q[qp+1],q[qp+2]));
        }
    }


    /** \brief Function for printing fatal error messages and exiting.
     *
     * Function for printing fatal error messages and exiting.
     * \param[in] p a pointer to the message to print.
     * \param[in] status the status code to return with. */
    public static void voro_fatal_error(String p,Config.Voropp status) {
        System.err.println("voro++: "+p);
        throw new RuntimeException("Status: "+status);
    }

/** \brief Prints a vector of positions.
 *
 * Prints a vector of positions as bracketed triplets.
 * \param[in] v the vector to print.
 * \param[in] fp the file stream to print to. */
    public static void voro_print_positions(DoubleArrayList v) {
        voro_print_positions(v, System.out);
    }
    public static void voro_print_positions(DoubleArrayList v, OutputStream fp) {
        if(v.size()>0) {
            try {
                fp.write(String.format("(%g,%g,%g)", v.getDouble(0), v.getDouble(1), v.getDouble(2)).getBytes());
                for (int k = 3; k < v.size(); k += 3) {
                    fp.write(String.format(" (%g,%g,%g)", v.getDouble(k), v.getDouble(k + 1), v.getDouble(k + 2)).getBytes());
                }
            }
            catch (IOException ex) {
                throw new RuntimeException(ex);
            }
        }
    }

    public static void voro_print_vector(IntArrayList v) {
        voro_print_vector(v, System.out);
    }
    public static void voro_print_vector(IntArrayList v,OutputStream fp) {
        int k,s=v.size();
        try {
            for (k = 0; k + 4 < s; k += 4) {
                fp.write(String.format("%g %g %g %g ", v.getInt(k), v.getInt(k + 1), v.getInt(k + 2), v.getInt(k + 3)).getBytes());
            }
            if (k + 3 <= s) {
                if (k + 4 == s)
                    fp.write(String.format("%g %g %g %g", v.getInt(k), v.getInt(k + 1), v.getInt(k + 2), v.getInt(k + 3)).getBytes());
                else fp.write(String.format("%g %g %g", v.getInt(k), v.getInt(k + 1), v.getInt(k + 2)).getBytes());
            } else {
                if (k + 2 == s) fp.write(String.format("%g %g", v.getInt(k), v.getInt(k + 1)).getBytes());
                else fp.write(String.format("%g", v.getInt(k)).getBytes());
            }
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }

    public static void voro_print_vector(DoubleArrayList v) {
        voro_print_vector(v, System.out);
    }
    public static void voro_print_vector(DoubleArrayList v,OutputStream fp) {
        int k,s=v.size();
        try {
            for (k = 0; k + 4 < s; k += 4) {
                fp.write(String.format("%g %g %g %g ", v.getDouble(k), v.getDouble(k + 1), v.getDouble(k + 2), v.getDouble(k + 3)).getBytes());
            }
            if (k + 3 <= s) {
                if (k + 4 == s)
                    fp.write(String.format("%g %g %g %g", v.getDouble(k), v.getDouble(k + 1), v.getDouble(k + 2), v.getDouble(k + 3)).getBytes());
                else fp.write(String.format("%g %g %g", v.getDouble(k), v.getDouble(k + 1), v.getDouble(k + 2)).getBytes());
            } else {
                if (k + 2 == s) fp.write(String.format("%g %g", v.getDouble(k), v.getDouble(k + 1)).getBytes());
                else fp.write(String.format("%g", v.getDouble(k)).getBytes());
            }
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }

    public static void voro_print_face_vertices(IntArrayList v) {
        voro_print_face_vertices(v, System.out);
    }
    public static void voro_print_face_vertices(IntArrayList v,OutputStream fp) {
        try {
            if (v.size() > 0) {
                int k = 0;
                int l = v.getInt(k);
                k++;
                if (l <= 1) {
                    if (l == 1) {
                        fp.write(String.format("(%d)", v.getInt(k)).getBytes());
                        k++;
                    } else {
                        fp.write("()".getBytes());
                    }
                } else {
                    int j = k + l;
                    fp.write(String.format("(%d", v.getInt(k)).getBytes());
                    k++;
                    for (; k < j; k++) fp.write(String.format(",%d", v.getInt(k)).getBytes());
                    fp.write(")\n".getBytes());
                }
                while (k < v.size()) {
                    l = v.getInt(k);
                    k++;
                    if (l <= 1) {
                        if (l == 1) {
                            fp.write(String.format(" (%d)", v.getInt(k)).getBytes());
                            k++;
                        } else fp.write(" ()\n".getBytes());
                    } else {
                        int j = k + l;
                        fp.write(String.format(" (%d", v.getInt(k)).getBytes());
                        k++;
                        for (; k < j; k++) fp.write(String.format(",%d", v.getInt(k)).getBytes());
                        fp.write(")\n".getBytes());
                    }
                }
            }
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }


}
