package simulate;
import java.io.*;

public class SpacePeriodicCubic extends Space {

    private transient double[] dr = new double[Space.D];
    protected final double[][] shift1 = new double[1][D]; //used by getOverflowShifts
    protected final double[][] shift3 = new double[3][D];
    
    
    public SpacePeriodicCubic() {
        super();
        periodic = true;
    }
    
    public final void uEr1Mr2(double[] u, double[] v1, double[] v2) {
        Space.uEv1Mv2(u, v1, v2);
        
        double dim = dimensions[0];
        double d2 = 0.5*dim;
        double uu = u[0];
        if(uu > 0) {u[0] = (uu > +d2) ? uu-dim : uu;}
        else       {u[0] = (uu < -d2) ? uu+dim : uu;}
        uu = u[1];
        if(uu > 0) {u[1] = (uu > +d2) ? uu-dim : uu;}
        else       {u[1] = (uu < -d2) ? uu+dim : uu;}
   /*     for(int i=Space.D; --i>=0; ) {  // i=0..D-1
           double dim = dimensions[i];
           if(u[i] > 0.5*dim) {
              u[i] -= dim;}
           else if(u[i] < -0.5*dim) {
              u[i] += dim;
           }
        }
    */
        return;
    }
    public final double r1Mr2_S(double[] v1, double[] v2) {
        uEr1Mr2(dr, v1, v2);
        return dr[0]*dr[0] + dr[1]*dr[1];
    }
    
    public final double r1iMr2i(int i, double[] v1, double[] v2) {
        double u = v1[i] - v2[i];
        double dim = dimensions[i];
        if(u > 0.5*dim) {
            u -= dim;}
        else if(u < -0.5*dim) {
            u += dim;
        }
        return u;
    }

    public void repositionMolecules() {
       for(AtomC a=(AtomC)parentPhase.firstAtom(); a!=null; a=(AtomC)a.getNextAtom()) {
           for(int i=Space.D; --i>=0; ) {
              if(a.r[i] < 0.0) {a.translate(i,dimensions[i]);}
              else if(a.r[i] > dimensions[i]) {a.translate(i,-dimensions[i]);}
           }
       }
    }
    
    
    /** Computes origins for periodic images
    */
    public double[][] imageOrigins(int nShells) {
        int nImages = (2*nShells+1)*(2*nShells+1)-1;
        double[][] origins = new double[nImages][D];
        int k = 0;
        for(int i=-nShells; i<=nShells; i++) {
            for(int j=-nShells; j<=nShells; j++) {
                if(i==0 && j==0) {continue;}
                origins[k][0] = i*dimensions[0];
                origins[k][1] = j*dimensions[1];
                k++;
            }
        }
        return origins;
    }

    /** Returns coordinate shifts needed to draw all images that overflow into central image
    */
    // 0, 1, or 3 shifts may be returned
    public double[][] getOverflowShifts(double[] r, double distance) {
        int shiftX = 0;
        int shiftY = 0;
        if(r[0]-distance < 0.0) {shiftX = +1;}
        else if(r[0]+distance > dimensions[0]) {shiftX = -1;}
        
        if(r[1]-distance < 0.0) {shiftY = +1;}
        else if(r[1]+distance > dimensions[1]) {shiftY = -1;}
        
        if(shiftX == 0) {
            if(shiftY == 0) {
                return shift0;
            }
            else {
                shift1[0][0] = 0.0;
                shift1[0][1] = shiftY*dimensions[1];
                return shift1;
            }
        }
        else { //shiftX != 0
            if(shiftY == 0) {
                shift1[0][0] = shiftX*dimensions[0];
                shift1[0][1] = 0.0;
                return shift1;
            }
            else {
                shift3[0][0] = shiftX*dimensions[0];
                shift3[0][1] = 0.0;
                shift3[1][0] = 0.0;
                shift3[1][1] = shiftY*dimensions[1];
                shift3[2][0] = shift3[0][0];
                shift3[2][1] = shift3[1][1];
                return shift3;
            }
        }
    }
}