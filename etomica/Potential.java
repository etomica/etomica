package etomica; 

/**
 * Basic methods needed to describe interaction between atoms.  
 * Molecular interactions are defined in the Potential1 and Potential2 classes by collecting
 * together Potential classes for all atoms pairs that can be made by taking one atom from each
 * molecule.
 *
 * @see Potential1
 * @see Potential2
 */
public abstract class Potential implements Simulation.Element, java.io.Serializable {
  
    public static String VERSION = "Potential:01.01.17";
    
    private final Simulation parentSimulation;
    private boolean added = false;
    private String name;
    
    public Potential(Simulation sim) {
        parentSimulation = sim;
    }
    
    public final Simulation parentSimulation() {return parentSimulation;}
    public final Class baseClass() {return Potential.class;}
    public final boolean wasAdded() {return added;}
    public final void setAdded(boolean b) {added = b;}
    public final String getName() {return name;}
    public final void setName(String name) {this.name = name;}
    
    /**
     * Returns the energy of interaction of the pair of atoms passed to the method
     */
    public abstract double energy(AtomPair pair);
          
    public double hyperVirial(AtomPair pair) {return 0.0;}
    /**
     * Returns true if the pair of atoms are considered by this potential to be overlapping in their current positions
     */
//    public abstract boolean overlap(AtomPair pair);

    /**
     * Returns the total (extensive) long-range correction to the energy, assuming g(r) = 1 beyond the truncation distance
     * Input are the number of atoms of each type (i.e., x1*N and x2*N, where x1 and x2
     * are the mole fractions of each species interacting according to this potential), and the volume
     * of the phase
     */
    public abstract double energyLRC(int n1, int n2, double V);
    
    //***************** end of methods for Potential interface *****************//
    
    /**
    * Methods needed to describe the behavior of a hard potential.  
    * A hard potential describes impulsive interactions, in which the energy undergoes a step
    * change at some separation (and perhaps orientation).  Atoms at this point are said to
    * collide.  Important examples of this type of potential are hard-sphere and square-well.
    *
    * @see Potential.Soft
    */
         
    public interface Hard {
            
    /**
    * Implements the collision dynamics.
    * The given pair of atoms are assumed to be at the point of collision.  This method is called
    * to change their momenta according to the action of the collision.  Extensions can be defined to
    * instead implement other, perhaps unphysical changes.
    */
        public void bump(AtomPair pair);
    /**
    * Computes the time of collision of the given pair, assuming no intervening collisions.
    * Usually assumes free-flight between collisions
    */ 
        public double collisionTime(AtomPair pair);
    /**
    * Returns the collision virial from the last collision processed by this potential
    * This quantity can be used to measure the pressure
    */
        public double lastCollisionVirial();
            
    /**
    *Returns the virial tensor from the last collision processed.  This is used to measure 
    *the pressure tensor, and eventually the surface tension
    */
        public etomica.Space.Tensor lastCollisionVirialTensor();
    }  //end of Potential.Hard

    /**
    * Methods needed to describe the behavior of a soft potential.  
    * A soft potential describes non-impulsive interactions, in which the energy at all points
    * has smooth, analytic behavior with no discontinuities.  
    *
    * @see Potential.Hard
    */
    public interface Soft {
        
        /**
        * Force exerted by one atom on the other.
        * By convention, this gives the change in the potential energy as atom1
        * of the pair is displaced as atom2 is held fixed, i.e., it is the
        * gradient of the potential with respect to the position of atom1
        *
        * @return the vector force that atom2 of the pair exerts on atom1
        */
        public Space.Vector force(AtomPair pair);
        /**
        * Returns the total long-range correction to the pressure
        *
        * @param n1 the number of atoms of species(Index)1, i.e., x1*N, where x1 is the mole fraction of species 1
        * @param n2 the number of atoms of species 2
        * @param V  the volume of the phase (needed to compute the density)
        */
        public double pressureLRC(int n1, int n2, double V);

        public double virial(AtomPair pair);
        public double hyperVirial(AtomPair pair);
        
        
    } //end of Potential.Soft
    
    public interface Reactive {
        
        public BondChangeData[] getBondChangeData();
        
        public static class BondChangeData {
            public Atom atom;
            public Atom[] oldPartners;
            public Atom[] newPartners;
            public Atom getAtom() {return atom;}
            public Atom[] getOldPartners() {return oldPartners;}
            public Atom[] getNewPartners() {return newPartners;}
        }
        
    }
        
}



