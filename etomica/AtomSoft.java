package simulate;

/**
 *  Atom suitable for handling impulsive collisions.
 *  Contains variables to keep track of collision partners and times.
 *  
 *  @author David Kofke
 *  @see Molecule
 */
public abstract class AtomSoft extends AtomC {

    /**
     * The components of the vector describing the force acting on the atom.
     */
    public final double[] f = new double[Space.D];
    
    /**
     * Flag indicating whether atom is subject to any forces
     * forceFree = true if no forces act on atom, and it moves in free flight
     * Default is true
     *
     * @see IntegratorHardField
     */
    private boolean forceFree;
    private final double[] partnerForce = new double[Space.D];

    public AtomSoft(Molecule parent, int index, double mass, double diameter) {
        super(parent, index, mass, diameter);
        this.zeroForce();
        setForceFree(true);
    }

    public final void setForceFree(boolean b) {forceFree = b;}
    public final boolean isForceFree() {return forceFree;}
    public final void zeroForce() {Space.uEa1(f, 0.0);}
    public final void addForce(double[] force) { Space.uPEv1(f, force);}
    public final void subtractForce(double[] force) {Space.uMEv1(f, force);}
    
  /*public final double[] getPartnerForce() {return partnerForce;}
  public final void setPartnerForce(double[] force) {
    Space.uEv1(partnerForce, force);
    this.time0 = 0.0;
  }
*/

}