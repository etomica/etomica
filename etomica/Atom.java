package simulate;
import java.awt.*;

/**
 *  Each instance of the class Atom holds the position and velocity of one
 *  physical atom; all simulation kinetics and dynamics are performed by 
 *  operating on these values.
 *  
 *  @author David Kofke
 *  @author C. Daniel Barnes
 *  @see Molecule
 */
public class Atom {

    /**
     * The Cartesian coordinates (Angstrom) describing the position of the atom.
     * Values range from zero to a maximum defined by the Space class.
     * @see Space#dimensions
     */
    public final double[] r = new double[Space.D];
    
    /**
     * The components of the momentum vector of the atom.
     * Units of momentum are (amu)(Angstrom)(ps)^-1
     */
    public final double[] p = new double[Space.D];
    
    /**
     * The components of the vector describing the force acting on the atom.
     * This field may later be placed in a planned <code>AtomSoft</code> 
     * subclass
     */
    public final double[] f = new double[Space.D];
    
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
    
    /**
     * Mass of atom in amu.
     * 1 amu = 1 atomic mass unit = (1/N_Avo) grams
     */
    double mass;
    
    /**
     * Reciprocal of the mass
     * @see #mass
     */
    double rm;
    
    /**
     * Color of the atom when drawn on the screen
     */
    Color color = Color.black;
    
    /**
     * Ratio of atom mass to mass of parent molecule.
     * Used to locate the center-of-mass of the molecule
     * @see Molecule#COM
     */
    double COMFraction;

    /**
     * Instance of molecule in which this atom resides.
     * Assigned in Atom constructor.
     * @see Molecule#makeAtoms
     */
    Molecule parentMolecule;
    
    /**
     * Identifier of atom within molecule.
     * Assigned by parent molecule when invoking Atom constructor.
     * @see Molecule#makeAtoms
     */
    int atomIndex;
    
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
    Atom collisionPartner;
    
    /**
     * next atom in linked list of atoms
     * @see #setNextAtom
     */
    Atom nextAtom;
    
    /**
     * previous atom in linked list of atoms
     *
     * @see #setNextAtom
     */
    Atom previousAtom;
    
    /**
     * (interatomic) potential governing interations between this atom and
     * its collision partner.
     */
    Potential collisionPotential;
    
    /**
     * Flag indicating whether atom is subject to any forces
     * forceFree = true if no forces act on atom, and it moves in free flight
     * Default is true
     *
     * @see IntegratorHardField
     */
    private boolean forceFree;
    
    private final double[] partnerForce = new double[Space.D];
    public final double[] rLast = new double[Space.D];
    
    /**
     * Constructs an atom with default values for mass, diameter, and color.
     * Default values are mass = 1.0 amu; diameter = 0.1 A; color = black.
     * Defaults for all coordinates and momenta are zero.
     *
     * @param parent       molecule in which atom resides
     * @param index        sequential index of atom as assigned by parent molecule
     */
    public Atom(Molecule parent, int index) {
        this(parent, index, 1.0, 0.1);
        /**   delete if everything works ok
        parentMolecule = parent;
        atomIndex = index;
        Space.uEa1(r,0.0);
        Space.uEa1(p,0.0);
        this.zeroForce();
        this.setMass(1.0);
        this.setDiameter(0.1);
        */
    }
    
    /**
     * Constructs an atom with passed values for mass, diameter and default color.
     * Default value of color is black.
     *
     * @param parent       molecule in which atom resides
     * @param index        sequential index of atom as assigned by parent molecule
     * @param mass         mass in amu
     * @param diameter     diameter (Angstrom) for drawing atom to screen
     */
    public Atom(Molecule parent, int index, double mass, double diameter) {
        parentMolecule = parent;
        atomIndex = index;
        Space.uEa1(r,0.0);
        Space.uEa1(p,0.0);
        this.zeroForce();
        this.setMass(mass);
        this.setDiameter(diameter);
        setForceFree(true);
    }
    
    public final Molecule getMolecule() {return parentMolecule;}
    
    public final int getSpeciesIndex() {return parentMolecule.getSpeciesIndex();}
    public final int getAtomIndex() {return atomIndex;}
    
    public final Color getColor() {return color;}
    public final void setColor(Color c) {this.color = c;}
    
    public final double getRm() {return rm;}
    
    public final void setForceFree(boolean b) {forceFree = b;}
    public final boolean isForceFree() {return forceFree;}
    /**
     * Sets reciprocal mass of this atom and updates mass accordingly.  Setting
     * rm to zero causes mass to be set to largest machine value.
     * 
     * @param rm   new value for reciprocal mass
     */
    public final void setRm(double rm) {
        this.rm = rm;
        mass = (rm==0.0) ? Double.MAX_VALUE : 1.0/rm;
    }
    
    public final double getMass() {return mass;}
    /**
     * Sets  mass of this atom and updates reciprocal mass accordingly.  Setting
     * mass to largest machine double (Double.MAX_VALUE) causes reciprocal mass 
     * to be set to zero.
     * 
     * @param mass   new value for mass
     */
    public final void setMass(double mass) {
        this.mass = mass;
        rm = (mass==Double.MAX_VALUE) ? 0.0 : 1.0/mass;
    }
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
    
    public final void zeroForce() {Space.uEa1(f, 0.0);}
    public final void addForce(double[] force) { Space.uPEv1(f, force);}
    public final void subtractForce(double[] force) {Space.uMEv1(f, force);}
    
//    public final void setR(double[] dr) {Space.uEv1(r,dr);}
//    public final void setR(int i, double dr) {r[i] = dr;}

    /**
     * Displaces atom by a vector dr.
     *
     * @param dr  vector specifying change in position
     */
    public final void translate(double[] dr) {Space.uPEv1(r,dr);}
    
    /**
     * Displaces atom the distance dr in one coordinate direction.
     *
     * @param i   index specifying coordinate direction of displacement
     * @param dr  displacement distance
     */
    public final void translate(int i, double dr) {r[i] += dr;}
    
    /**
     * Displaces atom by a vector dr, saving its original position,
     * which can be recovered by call to replace
     *
     * @param dr   vector specifying change in position
     * @see #replace
     */
     public final void displace(double[] dr) {
        Space.uEv1(rLast,r);
        translate(dr);
     }
     
     /**
      * Puts atom back in position held before last call to displace
      *
      * @see #displace
      */
      public final void replace() { Space.uEv1(r,rLast); }
    
     /**
      * Changes the momentum vector
      *
      * @param dp The change in the momentum vector
      */
    public final void accelerate(double[] dp) {Space.uPEv1(p,dp);}
 //   public final void accelerate(int i, double dp) {p[i] += dp;}
 //   public final void decelerate(double[] dp) {Space.uMEv1(p,dp);}
 //   public final void setP(double[] dp) {Space.uEv1(p,dp);}
 //   public final void setP(int i, double dp) {p[i] = dp;}
    
    public final void scaleP(int i, double scale) {p[i] *= scale;}
    
 
    /**
     * Computes and returns kinetic energy of atom.
     *
     * @return  kinetic energy in (amu)(Angstrom)<sup>2</sup>(ps)<sup>-2</sup>
     */
    public double kineticEnergy() {return 0.5*rm*Space.v1S(p);}
    
    /**
     * Computes and returns the potential energy of the atom due to its interactions
     * with all other atoms (intra and intermolecular) in the phase
     *
     * @return (potential energy)/kB in Kelvins
     */
    public final double potentialEnergy() {
        return intraPotentialEnergy() + interPotentialEnergy();
    }
    
    /**
     * Computes and returns the potential energy of the atom due to its interactions
     * with all other atoms in the molecule (intramolecular)
     *
     * @return (potential energy)/kB in Kelvins
     */
    public final double intraPotentialEnergy() {
        Potential1 p1 = parentMolecule.getP1();
        double energy = 0.0;
        for(Atom a=parentMolecule.firstAtom; a!=parentMolecule.lastAtom.getNextAtom(); a=a.getNextAtom()) {
            if(a == this) {continue;}
            energy += p1.getPotential(this,a).energy(this,a);
        }
        return energy;
    }
        
    /**
     * Computes and returns the potential energy of the atom due to its interactions
     * with all other atoms in all other molecules in the phase (intermolecular)
     *
     * @return (potential energy)/kB in Kelvins
     */
    public final double interPotentialEnergy() {
        Potential2[] p2 = parentMolecule.getP2();
        double energy = 0.0;
        Atom endAtom = parentMolecule.firstAtom;
        for(Atom a=parentMolecule.getPhase().firstAtom; a!=endAtom && a!=null; a=a.getNextAtom()) {
            energy += p2[a.getSpeciesIndex()].getPotential(this,a).energy(this,a);
            if(energy >= Double.MAX_VALUE) {return Double.MAX_VALUE;}
        }
        for(Atom a=parentMolecule.lastAtom.getNextAtom(); a!=null; a=a.getNextAtom()) {
            energy += p2[a.getSpeciesIndex()].getPotential(this,a).energy(this,a);
            if(energy >= Double.MAX_VALUE) {return Double.MAX_VALUE;}
        }
        return energy;
    }
    
    /**
     * Sets parameters relevant to next collision of atom.  Used in hard-potential
     * dynamics.
     *
     * @param time       time (ps) to next collision
     * @param partner    atom with which next collision will occur
     * @param potential  interatomic potential between these atoms
     */
    public final void setCollision(double time, Atom partner, Potential p) {
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
  
  public final Potential getCollisionPotential() {return collisionPotential;}
  
  public final Atom getCollisionPartner() {return collisionPartner;}
  
  /*public final double[] getPartnerForce() {return partnerForce;}
  public final void setPartnerForce(double[] force) {
    Space.uEv1(partnerForce, force);
    this.time0 = 0.0;
  }
*/
  public final Atom getNextAtom() {return nextAtom;}
  /**
   * Sets atom following this one in linked list, and sets this to be that
   * atom's previous atom in list
   * 
   * @param atom  the next atom in the list
   * @see Molecule#orderAtoms
   */
  public final void setNextAtom(Atom atom) {
    this.nextAtom = atom;
    if(atom != null) {atom.previousAtom = this;}
  }
  public final Atom getPreviousAtom() {return previousAtom;}
  
  /**
   * Computes and updates COMFraction as the ratio of this atom's mass to the
   * total mass of the parent molecule.
   */
  public final void updateCOMFraction() {COMFraction = mass/parentMolecule.getMass();}
  /**
   * Returns stored value of COMFraction.  Does not check to see if value needs
   * to be updated, but only returns stored value.
   *
   * @return   ratio of this atom's mass to the total mass of the parent molecule
   * @see #updateCOMFraction
   */
  public final double getCOMFraction() {return COMFraction;}

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
    
}