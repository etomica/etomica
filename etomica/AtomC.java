package simulate;

import java.awt.*;

/**
 *  Continuum Atom.  Atom in a continuous space (as opposed to a lattice).
 *  Holds the position and momentum of one physical atom; all simulation kinetics and dynamics are performed by 
 *  operating on these values.
 *  
 *  @author David Kofke
 *  @see Molecule
 */
public abstract class AtomC extends Atom {

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

    public final double[] rLast = new double[Space.D];
    
    /**
     * The components of the vector describing the force acting on the atom.
     */
    protected final double[] f = new double[Space.D];
    
    /**
     * Flag indicating whether atom is subject to any forces
     * forceFree = true if no forces act on atom, and it moves in free flight
     * Any call to method setForce sets forceFree to false; call to zeroForce sets it to true.
     * Default is true
     *
     * @see IntegratorHard
     */
    protected boolean forceFree;
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
     * Constructs an atom with passed values for mass, diameter and default color.
     * Default value of color is black.
     *
     * @param parent       molecule in which atom resides
     * @param index        sequential index of atom as assigned by parent molecule
     */
     
     AtomC nextAtomC;
     AtomC previousAtomC;
     
    public AtomC(Molecule parent, int index) {
        super(parent, index);
        if(parent != null) {    //null parent indicates atom is used only to generate other atoms
            Space.uEa1(r,0.0);
            Space.uEa1(p,0.0);
            this.setMass(1.0);
            setStationary(false);
            this.zeroForce();
        }
    
    }

    /**
     * Computes and returns the potential energy of the atom due to its interactions
     * with all other atoms in all other molecules in the phase (intermolecular)
     *
     * @return (potential energy)/kB in Kelvins
     */
    public double interPotentialEnergy() {
        Potential2[] p2 = parentMolecule.getP2();
        double energy = 0.0;
        Atom endAtom = parentMolecule.firstAtom;
        for(AtomC a=(AtomC)parentMolecule.getPhase().firstAtom(); a!=endAtom && a!=null; a=a.getNextAtomC()) {
            energy += p2[a.getSpeciesIndex()].getPotential(this,a).energy(this,a);
            if(energy >= Double.MAX_VALUE) {return Double.MAX_VALUE;}
        }
        for(AtomC a=(AtomC)parentMolecule.lastAtom.getNextAtom(); a!=null; a=a.getNextAtomC()) {
            energy += p2[a.getSpeciesIndex()].getPotential(this,a).energy(this,a);
            if(energy >= Double.MAX_VALUE) {return Double.MAX_VALUE;}
        }
        return energy;
    }
    
    /**
     * Computes and returns the potential energy of the atom due to its interactions
     * with all other atoms in the molecule (intramolecular)
     *
     * @return (potential energy)/kB in Kelvins
     */
    public final double intraPotentialEnergy() {
        if(parentMolecule.nAtoms == 1) return 0.0;
        Potential1 p1 = parentMolecule.getP1();
        double energy = 0.0;
        for(AtomC a=(AtomC)parentMolecule.firstAtom; a!=parentMolecule.lastAtom.getNextAtom(); a=a.getNextAtomC()) {
            if(a == this) {continue;}
            energy += p1.getPotential(this,a).energy(this,a);
        }
        return energy;
    }
        
    public final double getRm() {return rm;}
    
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
    
    public void setStationary(boolean b) {
        super.setStationary(b);
        if(b) Space.uEa1(p,0.0);   //zero momentum if set to stationary
    }

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
  
    public final void setForceFree(boolean b) {forceFree = b;}
    public final boolean isForceFree() {return forceFree;}
    public final void setForce(double[] force) {
        Space.uEv1(f, force);
        setForceFree(false);
    }
    public final double[] getForce() {return f;}
    public final void zeroForce() {
        Space.uEa1(f, 0.0);
        setForceFree(true);
    }
    public final void addForce(double[] force) { Space.uPEv1(f, force);}
    public final void subtractForce(double[] force) {Space.uMEv1(f, force);}
    
    public void setNextAtom(Atom atom) {
      super.setNextAtom(atom);
      this.nextAtomC = (AtomC)atom;
      if(atom != null) {((AtomC)atom).previousAtomC = this;}
    }
    
    public final AtomC getNextAtomC() {return nextAtomC;}
    public final AtomC getPreviousAtomC() {return previousAtomC;}
    
}