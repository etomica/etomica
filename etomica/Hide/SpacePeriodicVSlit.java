package simulate;
import java.io.*;

public class SpacePeriodicVSlit extends Space {

    private transient double[] dr = new double[Space.D];
    protected final double[][] shift1 = new double[1][D]; //used by getOverflowShifts
    
    
    public SpacePeriodicVSlit() {
        super();
        periodic = true;
    }
    
    public final void uEr1Mr2(double[] u, double[] v1, double[] v2) {
        Space.uEv1Mv2(u, v1, v2);
        int i = Space.D - 1;
           if(u[i] > 0.5*dimensions[i]) {
            u[i] -= dimensions[i];}
           else if(u[i] < -0.5*dimensions[i]) {
            u[i] += dimensions[i];}
        return;
    }
    public final double r1Mr2_S(double[] v1, double[] v2) {
        uEr1Mr2(dr, v1, v2);
        return dr[0]*dr[0] + dr[1]*dr[1];
    }

    public void repositionMolecules(Species species) {
       Atom nextSpeciesAtom = species.lastAtom.getNextAtom();
       int i = Space.D - 1;
       for (Atom a=species.firstAtom; a!=nextSpeciesAtom; a=a.getNextAtom()){
          if(a.r[i] < 0.0) {a.translate(i,dimensions[i]);}
          else if(a.r[i] > dimensions[i]) {a.translate(i,-dimensions[i]);}
       }
    }

  public void setScale(double s, int n) {   
      if(s>0 && n>=0) {
        scale = s/(double)(2*n+1);
        nShells = n;
        computeDrawSize();
      }
  }
    
    
    /** Computes origins for periodic images
    */
    protected void resetOrigins(int n) {
        nShells = n;
        nImages = 2*nShells;
        origins = new int[nImages][];
        int k = 0;
        for(int i=0; i<nImages; i++) {origins[i] = new int[D];}
        int i = 0;
        for(int j=-nShells; i<=nShells; i++) {
           if(j==0) {continue;}
           origins[k][0] = (int)(i*(drawSize[0]-1));
           origins[k][1] = (int)(j*(drawSize[1]-1));
           k++;
        }
    }

    /** Returns coordinate shifts needed to draw all images that overflow into central image
    */
    // 0, 1, or 3 shifts may be returned
    public double[][] getOverflowShifts(double[] r, double distance) {
        int shiftX = 0;
        int shiftY = 0;
        if(r[1]-distance < 0.0) {shiftY = +1;}
        else if(r[1]+distance > dimensions[1]) {shiftY = -1;}
        
        if(shiftY == 0) {return shift0;}
        else { //shiftY != 0
            shift1[0][1] = shiftY*dimensions[1];
            shift1[0][0] = 0.0;
            return shift1;
        }
    }
}