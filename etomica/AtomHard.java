package simulate;

/**
 *  Atom suitable for handling impulsive collisions.
 *  Contains variables to keep track of collision partners and times.
 *  
 *  @author David Kofke
 *  @see Molecule
 */
public abstract class AtomHard extends AtomC {

    /**
     * Time since last collision
     */
    public double time0;
    
    /**
     * Time (picoseconds) to next collision; used for hard-potential dynamics.
     * @see #collisionPartner
     * @see #setCollision
     */
    double collisionTime = Double.MAX_VALUE;
    
    /**
     * Atom that this atom will collide with next if no other collisions
     * change their trajectories before collision time; used for hard-potential
     * dynamics.
     * @see #collisionTime
     * @see #setCollision
     */
    AtomHard collisionPartner;

    /**
     * (interatomic) potential governing interations between this atom and
     * its collision partner.
     */
    PotentialHard collisionPotential;
    
    AtomHard nextAtomHard;
    AtomHard previousAtomHard;
    
    /**
     * Constructs an atom with no initialization if parent is null; otherwise constructs atom with default atomIndex = 0.  
     * Expected use of such an Atom is for the construction of other Atoms via makeAtom method
     */
    public AtomHard(Molecule parent) {
        this(parent,0);
    }
    /**
     * Constructs an atom with default values for mass, diameter, and color.
     * Default values are mass = 1.0 amu; diameter = 0.1 A.
     * Defaults for all coordinates and momenta are zero.
     *
     * @param parent       molecule in which atom resides
     * @param index        sequential index of atom as assigned by parent molecule
     */
    public AtomHard(Molecule parent, int index) {
        super(parent, index);
    }
    
    /**
     * Sets parameters relevant to next collision of atom.  Used in hard-potential
     * dynamics.
     *
     * @param time       time (ps) to next collision
     * @param partner    atom with which next collision will occur
     * @param potential  interatomic potential between these atoms
     */
    public final void setCollision(double time, AtomHard partner, PotentialHard p) {
        collisionTime = time;
        collisionPartner = partner;
        collisionPotential = p;
    }
    
  public final double getCollisionTime() {return collisionTime;}
  
  /**
   * Decreases time to next expected collision and increments time since 
   * last collision.
   *
   * @param interval   time interval
   * @see IntegratorHard#advanceAcrossTimeStep
   */
  public final void decrementCollisionTime(double interval) {
    this.collisionTime -= interval;
    this.time0 += interval;
  }
  
  public final PotentialHard getCollisionPotential() {return collisionPotential;}
  
  public final AtomHard getCollisionPartner() {return collisionPartner;}
  
  public void setNextAtom(Atom atom) {
    super.setNextAtom(atom);
    this.nextAtomHard = (AtomHard)atom;
    if(atom != null) {((AtomHard)atom).previousAtomHard = this;}
  }
    
  public final AtomHard getNextAtomHard() {return nextAtomHard;}
  public final AtomHard getPreviousAtomHard() {return previousAtomHard;}
}