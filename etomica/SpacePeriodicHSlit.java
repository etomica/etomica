package simulate;
import java.io.*;

public class SpacePeriodicHSlit extends Space {

    private transient double[] dr = new double[Space.D];
    protected final double[][] shift1 = new double[1][D]; //used by getOverflowShifts
    
    
    public SpacePeriodicHSlit() {
        super();
        periodic = true;
    }
    
    public final void uEr1Mr2(double[] u, double[] v1, double[] v2) {
        Space.uEv1Mv2(u, v1, v2);
        for(int i=Space.D-1; --i>=0; ) {  // i=0..D-2
           if(u[i] > 0.5*dimensions[i]) {
            u[i] -= dimensions[i];}
           else if(u[i] < -0.5*dimensions[i]) {
            u[i] += dimensions[i];}
        }
        return;
    }
    public final double r1Mr2_S(double[] v1, double[] v2) {
        uEr1Mr2(dr, v1, v2);
        return dr[0]*dr[0] + dr[1]*dr[1];
    }

    public void repositionMolecules() {
       for(Atom a=parentPhase.firstAtom(); a!=null; a=a.getNextAtom()) {
           for(int i=Space.D-1; --i>=0; ) {
              if(a.r[i] < 0.0) {a.translate(i,dimensions[i]);}
              else if(a.r[i] > dimensions[i]) {a.translate(i,-dimensions[i]);}
           }
       }
    }
    
    
    /** Computes origins for periodic images
    */
    public double[][] imageOrigins(int nShells) {
        int nImages = 2*nShells;
        double[][] origins = new double[nImages][D];
        int k = 0;
        int j = 0;
        for(int i=-nShells; i<=nShells; i++) {
           if(i==0) {continue;}
           origins[k][0] = i*(dimensions[0]-1);
           origins[k][1] = j*(dimensions[1]-1);
           k++;
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
        
        if(shiftX == 0) {return shift0;}
        else { //shiftX != 0
            shift1[0][0] = shiftX*dimensions[0];
            shift1[0][1] = 0.0;
            return shift1;
        }
    }
}