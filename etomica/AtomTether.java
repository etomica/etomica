package simulate;
import java.awt.*;

/**  
 *  @author David Kofke
 *  @see Atom
 */
public class AtomTether extends Atom {

    /**
     */
    public final double[] r0 = new double[Space.D];
    public final double[] r1 = new double[Space.D];
    public final double[] e = new double[Space.D];
    public final double[] f = new double[Space.D];
    public boolean bonded = false;
    public double stickingProbability;
    public double equilibriumLength;
    public double springConstant;
    
    public AtomTether(Molecule parent, int index, double[] r0, double length, double[] e) {
        super(parent, index);
        Space.uEv1(this.r0, r0);
        Space.uEv1Pa1Tv2(this.r1,r0,length,e);
    }
    
    public double[] force() {
        if(bonded) {
          Space.unitVector(e,r0,r1);
          double d = length() - equilibriumLength;
          Space.uEa1Tv1(f,springConstant*d,e);
        }
        else { Space.uEa1(f,0.0); }
        return f;
    }
    
    public double length() {return Space.v1Mv2_S(r0,r1);}
    
    
  /**
   * Draws this atom using current values of its position, diameter and color.
   * Drawing position is determined as follows.  The atoms coordinates in 
   * Angstroms are converted to pixels by applying a scaling factor; these
   * drawing coordinates may be shifted by some amount as given by the array
   * <code>origin</code> before the atom is drawn.
   *
   * @param g         graphics object to which atom is drawn
   * @param origin    origin of drawing coordinates (pixels)
   * @param scale     factor determining size of drawn image relative to
   *                  nominal drawing size
   */
  public void draw(Graphics g, int[] origin, double scale) {
    double toPixels = scale*Phase.TO_PIXELS;
    int sigmaP = (int)(toPixels*diameter);
    g.setColor(color);
    int x1 = origin[0] + (int)(toPixels*r0[0]);
    int y1 = origin[1] + (int)(toPixels*r0[1]);
    int x2 = origin[0] + (int)(toPixels*r1[0]);
    int y2 = origin[1] + (int)(toPixels*r1[1]);
    g.drawLine(x1,y1,x2,y2);
  }
    
}