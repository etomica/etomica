package simulate;
import java.io.*;
import java.awt.*;  //subclasses Component; for Graphics and Color in draw method
import java.beans.*;
import java.util.*;

/**
 * Each Phase must contain exactly one Space.  The dimensionality of the system 
 * is defined here.  The space object prescribes how
 * distances are measured, and it sets bounds on the allowable values of 
 * particle coordinates, defining what happens to molecules that go beyond
 * these bounds.  Thus subclasses of space can be defined to set up (for example)
 * periodic boundary conditions.
 *
 * Space also defines a large set of static utility methods to perform vector
 * arithmetic.  They are written specifically for vectors in a 2-dimensional
 * space.  Their purpose is to permit addition, multiplication, dot-products,
 * etc. of scalars and vectors without using a for-loop construct.  Extension
 * to simulation of 3-dimensional (or other) systems can be accomplished by
 * defining another Space routine with these methods appropriately redefined.
 * Since these are invoked as static methods, this object must replace (rather
 * than extend) the Space object defined here.  Methods used to compute vector
 * differences are written as instance (non-static) methods so they can
 * simply be overridden to re-define how distances are computed.  It is 
 * thought by us that the awkwardness of this approach is offset by the
 * speedup gained in using static, unrolled-loop methods to do these often-used
 * calculations.
 *
 * The utility methods are named according to a code that describes their function.
 *
 *      In the argument list, "u" represents the "left-hand side" if the 
 *      method calculation results in a vector; the array passed in this
 *      position will have its values updated according to the calculation.
 *      If the calculation results in a scalar, this is the return value
 *      of the method and no "u" is present in the argument list.
 *      
 *      v1, v2, etc. represent vectors appearing on the right-hand side
 *      of the vector equation, while a1, a2, etc. represent scalars 
 *      appearing there.
 *      
 *      Method names are coded to suggest how the vector equation
 *      would be written using standard notation.  Thus:
 *      
 *      uE --> "u equals"; returns a vector in the first argument
 *       P --> "plus"
 *       M --> "minus"
 *       T --> "times"
 *       D --> "dot"; inner product
 *       S --> "squared"; vector square magnitude
 *      PE --> "plus equals"  (+=); increments u
 *      ME --> "minus equals" (-=)
 *      underscores represent parentheses ("left paren" omitted if would be first character of method name)
 *
 * Space and its subclasses are Beans, and are directly accessed at design time.
 * At design time, a Space draws itself as a small yellow rectangle.  At run time it
 * may be invisible or it may be set to draw its boundaries.
 *
 * @see Phase
 */

public class Space extends Component {

 /**
  * Dimensionality of the space.  If this value is changed, all vector utility methods
  * must be updated accordingly.
  */
    public static final int D = 2;
 
 /**
  * True for spaces that have periodic boundaries
  */
    public boolean periodic;
 
 /**
  * Size of Phase (width, height) in Angstroms
  * Default value is 1.0 for each dimension.
  */
    protected final double[] dimensions = new double[D];
    
    private static final double[] work = new double[D];       //temporary work vector
 
 /**
  * @see SpacePeriodicCubic#getOverflowShifts
  */
    protected final double[][] shift0 = new double[0][D];
   
   /**
    * Phase in which this Space resides
    */
    protected Phase parentPhase;
   
   /**
    * Volume of the phase, in Angstroms^D.
    *
    * @see #computeVolume
    */
    public double volume;
      
    protected Random randomGenerator = new Random();
    
   /**
    * Creates Space with unit dimensions and computes its volume.
    */
    public Space() {
        dimensions[0] = 1.0;
        dimensions[1] = 1.0;
        computeVolume();
        periodic = false;
    }
    
    public Atom makeAtom(Molecule m, int index) {
        return new AtomCC(m, index);
    }
    
    private class AtomCC extends Atom {
        
        public Atom makeAtom(Molecule m, int i) {return null;}  //delete this
   
        public void draw(Graphics g, int[] origin, double scale) {}
        
        AtomCC(Molecule parent, int index) {
            super(parent, index);
            Space.uEa1(r,0.0);
            Space.uEa1(p,0.0);
            setMass(1.0);
            setStationary(false);
        }
        public final double[] r = new double[Space.D];  //Cartesian coordinates
        public final double[] p = new double[Space.D];  //Momentum vector
        private final double[] rLast = new double[Space.D];  //Displace/replace work vector

        /**
        * Flag indicating whether atom is subject to any forces
        * forceFree = true if no forces act on atom, and it moves in free flight
        * Any call to method setForce sets forceFree to false; call to zeroForce sets it to true.
        * Default is true
        *
        * @see IntegratorHard
        */
        protected boolean forceFree;
        double mass;  //Mass of atom in amu
        double rm;    //Reciprocal of the mass
         
        AtomCC nextAtomC;
        AtomCC previousAtomC;
         
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
            for(AtomCC a=(AtomCC)parentMolecule.firstAtom; a!=parentMolecule.lastAtom.getNextAtom(); a=a.getNextAtomCC()) {
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
          
        public void setNextAtom(Atom atom) {
          super.setNextAtom(atom);
          this.nextAtomC = (AtomCC)atom;
          if(atom != null) {((AtomCC)atom).previousAtomC = this;}
        }
            
        public final AtomCC getNextAtomCC() {return nextAtomC;}
        public final AtomCC getPreviousAtomCC() {return previousAtomC;}
    }
   /**
    * Sets parentPhase 
    *
    * @param p the phase to be identified as this Space's parentPhase
    * @see #setDimensions
    */
    public void setParentPhase(Phase p) {
        parentPhase = p;
    }
    
    public final boolean isPeriodic() {return periodic;}  
        
   /**
    * @return the dimensions[] array
    */
    public double[] getDimensions() {return dimensions;}
   
   /**
    * @return the value of dimensions[i]
    * @param i  index of desired value of dimensions[]
    */
    public double getDimensions(int i) {return dimensions[i];}
    
   /**
    * Sets one value in the dimensions array and updates drawSize array
    * and <code>volume</code>
    *
    * @param i   index of the value to be set
    * @param dim the new value of dimensions[i]
    */
    public void setDimensions(int i, double dim) {
        dimensions[i] = dim;
        computeVolume();
    }
   
   /**
    * Scales all dimensions by a constant multiplicative factor and recomputes volume
    *
    * @param scale the scaling factor. 
    */
    public void inflate(double scale) {
        Space.uTEa1(dimensions,scale);
        computeVolume();
    }
   
   /**
    * Computes the volume of the space based on the current values in
    * dimensions[], and assigns the result to <code>volume</code>.
    * Likely to override in subclasses.
    */
    protected void computeVolume() {
        volume = dimensions[0]*dimensions[1];
    }
        
   /** Set of vectors describing the displacements needed to translate the central image
    *  to all of the periodic images.  Returns a two dimensional array of doubles.  The
    *  first index specifies each perioidic image, while the second index indicates the
    *  x and y components of the translation vector.
    *  Likely to override in subclasses.
    *
    *  @param nShells the number of shells of images to be computed
    */
    public double[][] imageOrigins(int nShells) {
        return new double[0][D];
    }
       
   /**
    * @return shift0
    * @see #shift0
    * @see Phase#paint
    */
    public double[][] getOverflowShifts(double[] r, double distance) {return shift0;}  //called only if periodic
    
   /**
    * Design-time manipulation of dimensions[0].  Used in lieu of an array editor.
    */
    public void setXDimension(double dim) {setDimensions(0, dim);} //used for property list editor
    public double getXDimension() {return dimensions[0];}
    public void setYDimension(double dim) {setDimensions(1,dim);}
    public double getYDimension() {return dimensions[0];}
   
   /**
    * Method to handle placement of molecules that go beyond boundaries of Space.
    * Used by periodic-boundary subclasses.
    * 
    * @see Phase#paint
    */
    public void repositionMolecules() {return;}  

   /**
    * Method to draw Space at design time.  Draws a yellow rectangle to
    * confirm placement of space in phase.
    *
    * @param g graphic object to which image is drawn.
    */
    public void paint(Graphics g) {
        if(Beans.isDesignTime()) {
            g.setColor(Color.yellow);
            g.fillRect(0,0,getSize().width,getSize().height);
        }
    }
 
   /**
    * Draws a light gray outline of the space if <code>visible</code> is set to
    * <code>true</code>.  The size of the outline is determined by 
    * <code>drawSize[]</code>.
    *
    * @param g      the graphics object to which the image is drawn
    * @param origin the coordinate origin (in pixels) for drawing the image
    * @see Phase#paint
    * @see #computeDrawSize
    */
    public void drawFrame(Graphics g, int[] origin, double scale) {
        g.setColor(Color.gray.brighter());
        g.drawRect(origin[0],origin[1],(int)(scale*dimensions[0])-1,(int)(scale*dimensions[1])-1);
    }
    /* The remainder of the methods do vector arithmetic.
    
    */  
    
   /**
    * Assigns the constant a1 to all elements of the array u (u = a1).
    */
    public static void uEa1(double[] u, double a1) {
        u[0] = u[1] = a1;
    }
    
    public static void uEa1(int[] u, int a1) {
        u[0] = u[1] = a1;
    }
    
   /**
    * Element-by-element copy of array v1 to array u (u = v1)
    */
    public static void uEv1(double[] u, double[] v1) {
        u[0] = v1[0];
        u[1] = v1[1];
    }

   /**
    * Element-by-element copy of array v1 to array u (u = v1)
    */
    public static void uEv1(int[] u, int[] v1) {
        u[0] = v1[0];
        u[1] = v1[1];
    }
    
   /**
    * @return the scalar dot product of arrays v1 and v2
    */
    public static double v1Dv2(double[] v1, double[] v2) {
        return (v1[0]*v2[0] + v1[1]*v2[1]);
    }
   
   /**
    * @return the vector squared, the magnitude of the vector v1: v1^2 = v1 . v1
    */
    public static double v1S(double[] v1) {
        return (v1[0]*v1[0] + v1[1]*v1[1]);
    }
   
   /**
    * @return the scalar |v1-v2|^2
    */
    public static double v1Mv2_S(double[] v1, double[] v2) {
        double w0 = v1[0] - v2[0];
        double w1 = v1[1] - v2[1];
        return w0*w0 + w1*w1;
    }        
    
    public static void uEa1Tv1(double[] u, double a1, double[] v1) {
        u[0] = a1*v1[0];
        u[1] = a1*v1[1];
        return;
    }
    
    public static void uEa1Tv1(int[] u, double a1, double[] v1) {
        u[0] = (int)(a1*v1[0]);
        u[1] = (int)(a1*v1[1]);
        return;
    }
    
    public static void uEv1Mv2(double[] u, double[] v1, double[] v2) {
        u[0] = v1[0] - v2[0];
        u[1] = v1[1] - v2[1];
        return;
    }
    
    public void uEr1Mr2(double[] u, double[] v1, double[] v2) {  //for overriding by subclasses
        uEv1Mv2(u, v1, v2);
    }
    public double r1Mr2_S(double[] v1, double[] v2) {  //for overriding by subclasses
        return v1Mv2_S(v1, v2);
    }
    public double r1iMr2i(int i, double[] v1, double[] v2) {  //for overriding
        return v1[i] - v2[i];
    }
    
    public static void uPEv1(double[] u, double[] v1) {
        u[0] += v1[0];
        u[1] += v1[1];
        return;
    }
        
    public static void uMEv1(double[] u, double[] v1) {
        u[0] -= v1[0];
        u[1] -= v1[1];
        return;
    }
    
    public static void uPEa1(double[] u, double a1) {
        u[0] += a1;
        u[1] += a1;
    }
        
    public static void uMEa1(double[] u, double a1) {
        u[0] -= a1;
        u[1] -= a1;
    }
        
    public static void uTEa1(double[] u, double a1) {
        u[0] *= a1;
        u[1] *= a1;
        return;
    }
    
    public static void uDEa1(double[] u, double a1) {
        u[0] /= a1;
        u[1] /= a1;
    }
        
    public static void uPEa1Tv1(double[] u, double a1, double[] v1) {
        u[0] += a1*v1[0];
        u[1] += a1*v1[1];
        return;
    }
    
    public static void uMEa1Tv1(double[] u, double a1, double[] v1) {
        u[0] -= a1*v1[0];
        u[1] -= a1*v1[1];
        return;
    }
 
    public static void uEa1Tv1Ma2Tv2(double[] u, double a1, double[] v1, double a2, double[] v2) {
        u[0] = a1*v1[0] - a2*v2[0];
        u[1] = a1*v1[1] - a2*v2[1];
        return;
    }
    
    public static void uEa1Tv1Pa2Tv2(double[] u, double a1, double[] v1, double a2, double[] v2) {
        u[0] = a1*v1[0] + a2*v2[0];
        u[1] = a1*v1[1] + a2*v2[1];
        return;
    }
    
    public static void uEa1T_v1Mv2_(double[] u, double a1, double[] v1, double[] v2) {
        u[0] = a1*(v1[0]-v2[0]);
        u[1] = a1*(v1[1]-v2[1]);
    }
    
    public static void uEa1T_v1Mv2_(int[] u, double a1, int[] v1, int[] v2) {
        u[0] = (int)(a1*(v1[0]-v2[0]));
        u[1] = (int)(a1*(v1[1]-v2[1]));
    }
    
    public static void uEv1Pa1Tv2(double[] u, double[] v1, double a1, double[] v2) {
        u[0] = v1[0] + a1*v2[0];
        u[1] = v1[1] + a1*v2[1];
    }
    
    public static void uEv1Ma1Tv2(double[] u, double[] v1, double a1, double[] v2) {
        u[0] = v1[0] - a1*v2[0];
        u[1] = v1[1] - a1*v2[1];
    }
    
    public static void unitVector(double[] u, double[] v1, double[] v2) {
        uEv1Mv2(u, v1, v2);
        double norm = 1.0/v1S(u);
        uEa1Tv1(u,norm,u);
    }
    
    
    //would like to change name to randomPoint
    
    //random point in square centered on origin with each edge 2*rmax
    public static void randomVector(double[] u, double rmax, Random rand) {
        u[0] = (2.0*rand.nextDouble()-1.0)*rmax;
        u[1] = (2.0*rand.nextDouble()-1.0)*rmax;
    }
    
    public void randomVector(double[] u, Random rand) {  //random point in entire space
        u[0] = rand.nextDouble()*dimensions[0];
        u[1] = rand.nextDouble()*dimensions[1];
    }
    
    public double[] randomVector() {  //random point in entire space
        double[] u = new double[D];
        randomVector(u, randomGenerator);
        return u;
    }
}