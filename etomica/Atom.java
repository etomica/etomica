package etomica;
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
public class Atom implements Space.Occupant, java.io.Serializable {

    public static String getVersion() {return "01.01.17";}

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
        workVector = parentSimulation().space().makeVector();
        rLast = parentSimulation().space().makeVector();
        velocity = parentSimulation().space().makeVector();
        type = t; //this must precede makeCoordinate call
        coordinate = type.makeCoordinate(this);
      //  coordinate = Simulation.space.makeCoordinate(this);
        r = coordinate.position();
        p = coordinate.momentum();
        setStationary(false);
        if(atomLinkCount > 0) atomLink = new AtomLinker[atomLinkCount];
//        if(meterAgentCount > 0) meterAgents = new Meter.Agent(meterAgentCount);
    }
    public Atom(AtomGroup parent, AtomType t, int index) {
        parentGroup = parent;
        this.index = index;
        depth = 0;
        atomIndex = index;
        workVector = parentSimulation().space().makeVector();
        rLast = parentSimulation().space().makeVector();
        velocity = parentSimulation().space().makeVector();
        type = t; //this must precede makeCoordinate call
        coordinate = type.makeCoordinate(this);
      //  coordinate = Simulation.space.makeCoordinate(this);
        r = coordinate.position();
        p = coordinate.momentum();
        setStationary(false);
        if(atomLinkCount > 0) atomLink = new AtomLinker[atomLinkCount];
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
    
    public final AtomGroup parentGroup() {return parentGroup;}

    /**
     * Simulation in which this atom resides
     */
    public final Simulation parentSimulation() {return parentMolecule.parentSimulation();}
        
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
    public final void translateToRandom(etomica.Phase p) {translateTo(p.boundary().randomPosition());}
    public final void displaceBy(Space.Vector u) {rLast.E(r); translateBy(u);}
    public final void displaceBy(double d, Space.Vector u) {rLast.E(r); translateBy(d,u);}
    public final void displaceTo(Space.Vector u) {rLast.E(r); translateTo(u);}  
    public final void displaceWithin(double d) {workVector.setRandomCube(); displaceBy(d,workVector);}
    public final void displaceToRandom(etomica.Phase p) {rLast.E(r); translateToRandom(p);}
    public final void replace() {r.E(rLast);}
//    public final void inflate(double s) {r.TE(s);}

    public final void accelerateBy(Space.Vector u) {p.PE(u);}
    public final void accelerateBy(double d, Space.Vector u) {p.PEa1Tv1(d,u);}

    public final double kineticEnergy() {return coordinate.kineticEnergy();}
    public final void randomizeMomentum(double temperature) {  //not very sophisticated; random only in direction, not magnitude
        double magnitude = Math.sqrt(type.mass()*temperature*(double)parentSimulation().space().D());  //need to divide by sqrt(m) to get velocity
        p.setRandomSphere();
        p.TE(magnitude);
    }
    public final void scaleMomentum(double scale) {p.TE(scale);}

    public final Space.Vector position() {return r;}
    public final Space.Vector momentum() {return p;}
    public final double position(int i) {return r.component(i);}
    public final double momentum(int i) {return p.component(i);}
    public final Space.Vector velocity() {velocity.E(p); velocity.TE(type.rm()); return velocity;}  //returned vector is not thread-safe

    //needs work
    public final boolean preceeds(Atom atom) {
        if(atom == parentPhase().firstAtom()) return false;
        return true;
        /* //work this out later when atomGroup is in place
        //want to return true if atoms are the same atoms
        if(atom == null) return true;
        if(this.parentGroup == atom.parentGroup()) return this.index <= atom.index;//works also if both parentGroups are null
        else if(this.depth == atom.depth) return this.parentGroup.preceeds(atom.parentGroup());
        else if(this.depth < atom.depth) return this.preceeds(atom.parentGroup());
        else /*if(this.depth > atom.depth)* / return this.parentGroup.preceeds(atom);*/
    }
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
     
     public int index;
     public byte depth;
    
    private Atom nextAtom, previousAtom;
        
    public final Space.Coordinate coordinate;
    public final Space.Vector r, p;  //position, momentum
    private final Space.Vector workVector;
    private final Space.Vector rLast;
        
    public AtomType type;
    private final Space.Vector velocity;
    private AtomGroup parentGroup;
    
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
    public AtomLinker[] atomLink;
    
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

    public static int count(AtomIterator iter) {
        int n = 0;
        iter.reset();
        while(iter.hasNext()) {iter.next(); n++;}
        return n;
    }
 
}