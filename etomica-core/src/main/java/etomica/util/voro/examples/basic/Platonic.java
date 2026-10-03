// Platonic solids example code
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro.examples.basic;

import etomica.util.voro.VoronoiCell;

public class Platonic {

    public static final double Phi=0.5*(1+Math.sqrt(5.0));
    public static final double phi=0.5*(1-Math.sqrt(5.0));


    public static void main(String[] args) {
        VoronoiCell v = new VoronoiCell();

        // Create a tetrahedron
        v.init(-2,2,-2,2,-2,2);
        v.plane(1,1,1);
        v.plane(1,-1,-1);
        v.plane(-1,1,-1);
        v.plane(-1,-1,1);

        v.draw_gnuplot(0,0,0,"tetrahedron.gnu");

        // Create a cube. Since this is the default shape
        // we don't need to do any plane cutting.
        v.init(-1,1,-1,1,-1,1);
        v.draw_gnuplot(0,0,0,"cube.gnu");

        // Create an octahedron
        v.init(-2,2,-2,2,-2,2);
        v.plane(1,1,1);
        v.plane(-1,1,1);
        v.plane(1,-1,1);
        v.plane(-1,-1,1);
        v.plane(1,1,-1);
        v.plane(-1,1,-1);
        v.plane(1,-1,-1);
        v.plane(-1,-1,-1);

        v.draw_gnuplot(0,0,0,"octahedron.gnu");

        // Create a dodecahedron
        v.init(-2,2,-2,2,-2,2);
        v.plane(0,Phi,1);
        v.plane(0,-Phi,1);
        v.plane(0,Phi,-1);
        v.plane(0,-Phi,-1);
        v.plane(1,0,Phi);
        v.plane(-1,0,Phi);
        v.plane(1,0,-Phi);
        v.plane(-1,0,-Phi);
        v.plane(Phi,1,0);
        v.plane(-Phi,1,0);
        v.plane(Phi,-1,0);
        v.plane(-Phi,-1,0);

        v.draw_gnuplot(0,0,0,"dodecahedron.gnu");

        // Create an icosahedron
        v.init(-2,2,-2,2,-2,2);
        v.plane(1,1,1);
        v.plane(-1,1,1);
        v.plane(1,-1,1);
        v.plane(-1,-1,1);
        v.plane(1,1,-1);
        v.plane(-1,1,-1);
        v.plane(1,-1,-1);
        v.plane(-1,-1,-1);
        v.plane(0,phi,Phi);
        v.plane(0,phi,-Phi);
        v.plane(0,-phi,Phi);
        v.plane(0,-phi,-Phi);
        v.plane(Phi,0,phi);
        v.plane(Phi,0,-phi);
        v.plane(-Phi,0,phi);
        v.plane(-Phi,0,-phi);
        v.plane(phi,Phi,0);
        v.plane(phi,-Phi,0);
        v.plane(-phi,Phi,0);
        v.plane(-phi,-Phi,0);

        v.draw_gnuplot(0,0,0,"icosahedron.gnu");

    }
}
