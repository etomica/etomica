package simulate;

/**
 *  Atom suitable for handling impulsive collisions.
 *  Contains variables to keep track of collision partners and times.
 *  
 *  @author David Kofke
 *  @see Molecule
 */
public abstract class AtomSoft extends AtomC {

    private final double[] partnerForce = new double[Space.D];

    public AtomSoft(Molecule parent, int index) {
        super(parent, index);
    }

  /*public final double[] getPartnerForce() {return partnerForce;}
  public final void setPartnerForce(double[] force) {
    Space.uEv1(partnerForce, force);
    this.time0 = 0.0;
  }
*/

}