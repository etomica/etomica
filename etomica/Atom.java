package simulate;
import java.awt.*;
import java.awt.event.ActionEvent;

 /**
  * Each instance of the class Atom holds the coordinates of one physical atom; 
  * all simulation kinetics and dynamics are performed by operating on these values.
  * The coordinate is provided the the Space object when the atom is instantiated.
  * In addition, an Atom has a type that holds information such as how it is drawn,
  * and possibly information to be used as parameters in a potential.  Each atom
  * holds an Integrator.Agent object that may be used to store information needed
  * by the integrator.
  * Atoms are constructed when Species creates a Molecule.
  * 
  * @author David Kofke
  * @author C. Daniel Barnes
  * @see Molecule
  */
public final class Atom implements Space.Occupant, java.io.Serializable {

    /**
     * Constructs an atom with default values for mass, diameter, and color.
     * Defaults for all coordinates and momenta are zero.
     * 
     * @param parent molecule in which atom resides
     * @param index sequential index of atom as assigned by parent molecule
     * @param t the type of the atom
     */
    public Atom(Molecule parent, AtomType t, int index) {
        parentMolecule = parent;
        atomIndex = index;
        workVector = parentPhase().parentSimulation().space().makeVector();
        rLast = parentPhase().parentSimulation().space().makeVector();
        velocity = parentPhase().parentSimulation().space().makeVector();
        type = t; //this must precede makeCoordinate call
        coordinate = type.makeCoordinate(this);
      //  coordinate = Simulation.space.makeCoordinate(this);
        r = coordinate.position();
        p = coordinate.momentum();
        setStationary(false);
        if(atomLinkCount > 0) atomLink = new Atom.Linker[atomLinkCount];
//        if(meterAgentCount > 0) meterAgents = new Meter.Agent(meterAgentCount);
    }
                    
    /**
     * Assigns the atom's integrator agent to the given instance
     * 
     * @param ia
     */
    public void setIntegratorAgent(Integrator.Agent ia) {this.ia = ia;}
        
    /**
     * Molecule that contains this atom
     */
    public final Molecule parentMolecule() {return parentMolecule;}

    /**
     * Phase in which this atom resides
     */
    public final Phase parentPhase() {return parentMolecule.parentPhase();}
        
    /**
     * Coordinates of this atom.
     * When the atom is constructed the coordinate class is provided by the 
     * governing Space for the simulation.
     */
    public final Space.Coordinate coordinate() {return coordinate;}
            
    /**
     * Index identifying the species of this atom's molecule.
     */
    public final int speciesIndex() {return parentMolecule.speciesIndex();}

    /**
     * Integer assigned to this atom by its parent molecule.
     * Assigned during construction of atom.
     */
    public final int atomIndex() {return atomIndex;}
        
    /**
     * The color of the atom when drawn to the screen.
     * The atom color is often decided by the ColorScheme class, via the atom's setColor method.
     * If the ColorScheme does nothing, the color is given by the atom's type field.
     * 
     * @see ColorScheme
     */
    public final Color getColor() {return color;}

    /**
     * Sets the color of the atom for drawing to the screen.
     * 
     * @param c
     * @see ColorScheme
     * @see getColor
     */
    public final void setColor(Color c) {this.color = c;}

        
    /**
     * Sets the atom to be stationary or movable.
     * The atom does not enforce the condition of being stationary, in that it does not
     * override attempts to move it if it is set as stationary.  Rather, this is a flag
     * that can be set and checked by an integrator when deciding how or whether to 
     * move the atom.
     * 
     * @param b If true, the atom becomes stationary, if false it can be moved.
     */
    public void setStationary(boolean b) {
        stationary = b;
        if(!stationary) scaleMomentum(0.0);
    }

    /**
     * Flag to indicate of the atom can or cannot be moved
     * 
     * @see setStationary
     */
    public final boolean isStationary() {return stationary;}

    /**
     * Sets atom following this one in linked list, and sets this to be that atom's previous atom in list
     * 
     * @param atom the next atom in the list
     */
    public final void setNextAtom(Atom atom) {
        nextAtom = atom;
        if(atom!=null) {atom.previousAtom = this;}
    }

    /**
     * Sets this atom's previous atom to be null
     * 
     * @see setNextAtom
     */
    public final void clearPreviousAtom() {previousAtom = null;}

    /**
     * @return the next atom in the linked list of atoms
     * @see setNextAtom
     */
    public final Atom nextAtom() {return nextAtom;}

    /**
     * @return the previous atom in the linked list of atoms
     * @see setNextAtom
     */
    public final Atom previousAtom() {return previousAtom;} 
    
    /**
     * @return mass of the atom, in Daltons
     */
    public final double mass() {return type.mass();}

    /**
     * @return reciprocal of the mass of the atom
     */
    public final double rm() {return type.rm();}

    /**
     * Draws the atom to the given graphics object
     * 
     * @param g
     * @param origin
     * @param scale
     */
    public void draw(Graphics g, int[] origin, double toPixels) {type.draw(g, origin, toPixels, this);}

    /**
     * Moves the atom by some vector distance
     * 
     * @param u
     */
    public final void translateBy(Space.Vector u) {r.PE(u);}
    /**
     * Moves the atom by some vector distance
     * 
     * @param u
     */
    public final void translateBy(double d, Space.Vector u) {r.PEa1Tv1(d,u);}
    /**
     * Moves the atom by some vector distance
     * 
     * @param u
     */
    public final void translateTo(Space.Vector u) {r.E(u);}      
    public final void translateToRandom(simulate.Phase p) {translateTo(p.boundary().randomPosition());}
    public final void displaceBy(Space.Vector u) {rLast.E(r); translateBy(u);}
    public final void displaceBy(double d, Space.Vector u) {rLast.E(r); translateBy(d,u);}
    public final void displaceTo(Space.Vector u) {rLast.E(r); translateTo(u);}  
    public final void displaceWithin(double d) {workVector.setRandomCube(); displaceBy(d,workVector);}
    public final void displaceToRandom(simulate.Phase p) {rLast.E(r); translateToRandom(p);}
    public final void replace() {r.E(rLast);}
//    public final void inflate(double s) {r.TE(s);}

    public final void accelerateBy(Space.Vector u) {p.PE(u);}
    public final void accelerateBy(double d, Space.Vector u) {p.PEa1Tv1(d,u);}

    public final double kineticEnergy() {return coordinate.kineticEnergy();}
    public final void randomizeMomentum(double temperature) {  //not very sophisticated; random only in direction, not magnitude
        double magnitude = Math.sqrt(type.mass()*temperature*(double)parentPhase().parentSimulation().space().D());  //need to divide by sqrt(m) to get velocity
        p.setRandomSphere();
        p.TE(magnitude);
    }
    public final void scaleMomentum(double scale) {p.TE(scale);}

    public final Space.Vector position() {return r;}
    public final Space.Vector momentum() {return p;}
    public final double position(int i) {return r.component(i);}
    public final double momentum(int i) {return p.component(i);}
    public final Space.Vector velocity() {velocity.E(p); velocity.TE(type.rm()); return velocity;}  //returned vector is not thread-safe

    public Integrator.Agent ia;
        
    /**
    * Color of the atom when drawn on the screen
    * This color is set by the colorScheme object in the atom's species
    */
    Color color = Default.ATOM_COLOR;
        
    /**
    * Flag indicating whether atom is stationary or mobile.
    * Default is false (atom is mobile)
    */
    private boolean stationary;
       
    /**
    * Instance of molecule in which this atom resides.
    * Assigned in Atom constructor.
    */
//    final Molecule parentMolecule;
     Molecule parentMolecule;
        
    /**
    * Identifier of atom within molecule.
    * Assigned by parent molecule when invoking Atom constructor.
    */
//    final int atomIndex;
     int atomIndex;
    
    private Atom nextAtom, previousAtom;
        
    public final Space.Coordinate coordinate;
    public final Space.Vector r, p;  //position, momentum
    private final Space.Vector workVector;
    private final Space.Vector rLast;
        
    public AtomType type;
    private final Space.Vector velocity;
    
 /** This is an array of AtomLinkers that enable lists of atoms to be constructed
  *  and associated with this atom.  Useful for setting up neighbor lists, for example, or
  *  for defining bonds between atoms.
  *  Each element of this array points to the first linker in some linked list of atoms.
  *  Different simulation elements may wish to establish lists; they can do
  *  so without interfering with each other by requesting an index from the requestAtomLinkIndex()
  *  method of this class.  This returns a unique index that the element can use to construct a
  *  linked list of atoms beginning from atomLink[index].
  *  In many situations, this array is not used at all.
  *  @see Iterator.FixedNeighbors
  */
    public Atom.Linker[] atomLink;
    
    /**
     * Counter of the number of atom link index requests that have been fielded so far.
     */
    private static int atomLinkCount = 0;
    
    /**
     * Returns a unique index that can be used by a simulation element to set up a linked
     * list of atoms associated with each atom.  If idx is the value returned by this
     * method, the element is permitted to set up a linked list of atoms beginning
     * with atomLink[idx].
     */
    public static int requestAtomLinkIndex() {atomLinkCount++; return atomLinkCount-1;}
    
//    public Meter.Agent[] meterAgents;
    private static int meterAgentCount = 0;
    public static int requestMeterAgentIndex() {
        meterAgentCount++; 
        return meterAgentCount-1;
    }
 
 //*********************** ATOM ITERATORS *************************//
 
 /**
  * Base interface for atom iterators.
  * Atom iterators yield a sequence of atoms in successive calls to the next() method.
  * When all atoms have been returned, hasNext() returns false.
  * Iterators are often defined to progress "Up" or "Down" the set of atoms.
  * "Up" and "Down" are arbitrary designations, except that the iterators guarantee
  * that if atom 1 is "up list" of atom 2, then atom 2 is "down list" of atom 1.
  * "Up" and "Down" relation between any atoms may change during the course of the 
  * simulation, but at any instant the order is consistent and reproducible.
  * "Neighbor" iterators yield only atoms that are considered to be "neighbors" of
  * a specified atom.  The definition of "neighbor" depends on the iterator.  "Up neighbors"
  * are those neighbors uplist of the atom; likewise with "Down neighbors".
  *
  * @see IteratorFactory
  */
    public interface Iterator extends java.io.Serializable {        

        /**
         * @return true if the iterator will return another atom with a subsequent call to next(), false otherwise
         */
        public boolean hasNext();

        /**
         * @return the next atom in the list
         */
        public Atom next();

        /**
         * Resets the iterator in reference to the given atom.
         * Exactly how the given atom affects the reset depends on the particular iterator.
         * 
         * @param a
         */
        public void reset(Atom a);

        /**
         * Resets the iterator, so that it is ready to go through its list again.
         */
        public void reset();

        /**
         * Performs the given Action on each atom in the list in sequence.
         * 
         * @param act
         * @see Atom.Action
         */
        public void allAtoms(AtomAction act);
            
//        public Phase phase();
//        public Atom first();  //simply returns the first atom, without affecting status of iterator

    /**
     * Generic iterator that permits addition and removal of atoms.
     */
       public static final class List implements Iterator {
           private Atom.Linker first, next;
           private boolean hasNext;
           public List() {hasNext = false;}
           public boolean hasNext() {return hasNext;}
           public void reset() {
              next = first;
              hasNext = (next != null);
           }
           /**
            * reset with atom argument not properly implemented yet
            */
           public void reset(Atom a) {reset();}
           public Atom next() { //does not check that next is non-null
              Atom atom = next.atom();
              next = next.next();
              hasNext = (next != null);
              return atom;
           }
           public void allAtoms(AtomAction act) {
              for(Atom.Linker link=first; link!=null; link=link.next()) {
                  act.actionPerformed(link.atom());
              }
           }
           /**
            * Adds an atom to the set of atoms given by this iterator
            */
           public void addAtom(Atom a) {
              first = new Atom.Linker(a,first);
           }
           //will someday add a removeAtom method
           
       }//end of Atom.Iterator.List
     
    /**
     * Iterator that expires after returning a single atom
     */
        public static final class Singlet implements Iterator {
            private Atom atom;
            private boolean hasNext;
            public Singlet() {hasNext = false;}
            public Singlet(Atom a) {reset(a);}
            public boolean hasNext() {return hasNext;}
            public void reset() {hasNext = (atom != null);}
            public void reset(Atom a) {atom = a; reset();}
            public Atom next() {hasNext = false; return atom;}
            public void allAtoms(AtomAction act) {act.actionPerformed(atom);}
        }
        
        /**
         * Iterator that progresses up a list of atoms.
         * Order of atoms is that given by the linked list of atoms, which changes
         * only if atoms are added or removed to/from the phase.
         */
        public static class Up implements Iterator {
            protected Atom atom, nextAtom;
            protected boolean hasNext;
            private Phase phase;
            /**
             * Sets phase but does not reset iterator.  
             * Initial state is hasNext = false.
             */
            public Up(Phase p) {phase = p; hasNext = false;}
            /**
             * Sets phase and resets iterator using given atom.  
             * Initial state is hasNext = true if atom is not null.
             */
            public Up(Phase p, Atom a) {phase = p; reset(a);}
            public Phase phase() {return phase;}
            public boolean hasNext() {return hasNext;}
            /**
             * Sets the iterator so the next atom is the one given (which may be null)
             */
            public void reset(Atom a) {
                atom = a;
                hasNext = (a != null);
            }
            /**
             * Sets the iterator so that the next atom is the first atom of the phase
             */
            public void reset() {reset(phase.firstAtom());}
            public Atom next() {
                nextAtom = atom;
                atom = atom.nextAtom();
                hasNext = (atom != null);
                return nextAtom;
            }
            /**
             * Performs the given action on all atoms in the phase.  Unaffected by any prior calls to reset.
             */
            public void allAtoms(AtomAction act) {
                for(Atom a=phase.firstAtom(); a!=null; a=a.nextAtom()) {act.actionPerformed(a);}
            }
        } //end of Atom.Iterator.Up
        
        /**
         * Iterates over all neighbors uplist from given atom.
         * This iterator defines <i>all</i> atoms in phase as neighbors, so it
         * simply extends the Up iterator, modifying the reset method to start with
         * the atom following the given atom, rather than the given atom itself.
         * Also, the no-argument reset performs a reset using the current atom, rather than 
         * resetting to neighbor of first atom in phase.
         */
        public static final class UpNeighbor extends Up {
            private Atom first;
            public UpNeighbor(Phase p) {super(p);}
            public UpNeighbor(Phase p, Atom a) {super(p,a);} 
            /**
             * Resets iterator so that the next atom is the one just upList of the given atom.
             */
            public void reset(Atom a) {
                atom = a;
                if(a == null) {hasNext = false; return;}
                first = a.nextAtom();
                super.reset(first);
            }
            /**
             * Resets iterator to the condition it was in after the last call to reset(Atom a).
             * This will be hasNext = false if reset(Atom a) was not called previously.
             */
            public void reset() {super.reset(first);}
            /**
             * Performs the given action on all atoms uplist of the one indicated in the last call to reset(Atom).  
             * If reset has not been called before, performs no action.
             */
            public void allAtoms(AtomAction act) {
                if(first == null) return;
                for(Atom a=first; a!=null; a=a.nextAtom()) {act.actionPerformed(a);}
            }
        }
        
        /**
         * Iterator that progresses down a list of atoms.
         * Order of atoms is that given by the linked list of atoms, which changes
         * only if atoms are added or removed from the phase.
         */
        public static class Down implements Iterator {
            protected Atom atom;
            protected boolean hasNext;
            private Phase phase;
            /**
             * Sets phase but does not reset iterator.  
             * Initial state is hasNext = false.
             */
            public Down(Phase p) {phase = p; hasNext = false;}
            /**
             * Sets phase and resets iterator using given atom.  
             * Initial state is hasNext = true if atom is not null.
             */
            public Down(Phase p, Atom a) {phase = p; reset(a);}
            public Phase phase() {return phase;}
            public boolean hasNext() {return hasNext;}
            public void reset(Atom a) {
                atom = a;
                hasNext = (a != null);
            }
            /**
             * Resets iterator to the first atom of the list.
             * Iterator will return only this atom and then expire, since there is nothing downlist of it
             */
            public void reset() {reset(phase.firstAtom());}
            public Atom next() {
                Atom nextAtom = atom;
                atom = atom.previousAtom();
                if(atom == null) {hasNext = false;}
                return nextAtom;
            }
            /**
             * Performs the given action on all atoms in the phase, starting from the last to the first.
             */
            public void allAtoms(AtomAction act) {
                for(Atom a=phase.lastAtom(); a!=null; a=a.previousAtom()) {act.actionPerformed(a);}
            }
        } //end of Atom.Iterator.Down
        
        /**
         * Iterates over all neighbors downlist from given atom.
         * This iterator defines <i>all</i> atoms in phase as neighbors, so it
         * simply extends the Down iterator, modifying the reset method to start with
         * the atom preceding the given atom, rather than the given atom itself.
         * Also, the no-argument reset performs a reset using the current atom, rather than 
         * resetting to neighbor of first atom in phase.
         */
        public static final class DownNeighbor extends Down {
            private Atom first;
            public DownNeighbor(Phase p) {super(p);}
            public DownNeighbor(Phase p, Atom a) {super(p,a);}
            /**
             * Resets iterator so that the next atom is the one just downList of the given atom
             */
            public void reset(Atom a) {
                atom = a;
                if(a == null) {hasNext = false; return;}
                first = a.previousAtom();
                super.reset(first);
            }
            /**
             * Resets iterator to the condition it was in after the last call to reset(Atom a)
             * This will be hasNext = false if reset(Atom a) was not called previously.
             */
            public void reset() {super.reset(first);}
            /**
             * Performs the given action on all atoms uplist of the one indicated in the last call to reset(Atom).  
             * If reset has not been called before, performs no action.
             */
            public void allAtoms(AtomAction act) {
                if(first == null) return;
                for(Atom a=first; a!=null; a=a.previousAtom()) {act.actionPerformed(a);}
            }
        }
    }
    
    /**
     * Class for constructing linked lists of Atoms.
     * Each Linker points to one atom and another Linker, the next one in the list.
     * Although each atom has built-in ability to link to one next and one previous atom, these
     * Linkers are needed to construct other lists of atoms, particularly for neighbor lists.
     */
    public static class Linker implements java.io.Serializable {
        private final Atom atom;
        private Atom.Linker next = null;
        //Constructors
        public Linker(Atom a) {atom = a;}
        public Linker(Atom a, Linker l) {atom = a; next = l;}
        //Access methods
        public final Atom atom() {return atom;}
        public final Atom.Linker next() {return next;}
        public final void setNext(Atom.Linker l) {next = l;}
    }
}