package simulate;
import java.io.*;
import java.awt.*;//for Graphics
import java.util.*; //for neighborList
import java.beans.Beans;

public class AtomHardDisk extends AtomHard implements AtomDisk {

    /**
     * Diameter (Angstrom) used to specify the size of the atom when drawn to screen.
     * Not necessarily related to collision diameter in interatomic potential.
     * @see PotentialHardDisk#collisionDiameter
     */
    double diameter;
    
    /**
     * One-half of diameter
     * @see #diameter
     */
    double radius;
    
    AtomHardDisk nextAtomHardDisk;
    AtomHardDisk previousAtomHardDisk;
    
    /**
     * Constructs an atom with no initialization if parent is null; otherwise constructs atom with default atomIndex = 0.  
     * Expected use of such an Atom is for the construction of other Atoms via makeAtom method
     */
    public AtomHardDisk(Molecule parent) {
        this(parent,0);
    }
    
    public AtomHardDisk(Molecule parent, int index) {
        super(parent, index);
        if(parent != null) this.setDiameter(0.1);
    }
    
    public Atom makeAtom(Molecule m, int i) {return new AtomHardDisk(m,i);}

    public final double getDiameter() {return diameter;}
    /**
     * Sets diameter of this atom and updates radius accordingly.
     *
     * @param d   new value for diameter
     */
    public final void setDiameter(double d) {
        diameter = d;
        radius = 0.5*d;
    }
        
    public final double getRadius() {return radius;}
    /**
     * Sets radius of this atom and updates diameter accordingly.
     *
     * @param r   new value for radius
     */
    public final void setRadius(double r) {this.setDiameter(2.0*r);}
    
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
//    parentMolecule.parentSpecies.colorScheme.setAtomColor(this);
    g.setColor(color);
    int xP = origin[0] + (int)(toPixels*(r[0]-radius));
    int yP = origin[1] + (int)(toPixels*(r[1]-radius));
    g.fillOval(xP,yP,sigmaP,sigmaP);
  }
  
  
  public void setNextAtom(Atom atom) {
    super.setNextAtom(atom);
    this.nextAtomHardDisk = (AtomHardDisk)atom;
    if(atom != null) {((AtomHardDisk)atom).previousAtomHardDisk = this;}
  }
    
  public final AtomHardDisk getNextAtomHardDisk() {return nextAtomHardDisk;}
  public final AtomHardDisk getPreviousAtomHardDisk() {return previousAtomHardDisk;}
}