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
  */
public class Atom implements java.io.Serializable {

    public static String getVersion() {return "Atom:01.07.12";}
    
    public Atom(AtomGroup parent, AtomType t) {
        parentGroup = parent;
        type = t;
        coord = parentSimulation().space().makeCoordinate(this);//must follow setting of type field
        if(atomLinkCount > 0) atomLink = new AtomLinker[atomLinkCount];
    }
                    
    /**
     * Assigns the atom's integrator agent to the given instance
     * 
     * @param ia
     */
    public void setIntegratorAgent(Integrator.Agent ia) {this.ia = ia;}
            
    public final AtomGroup parentGroup() {return parentGroup;}
    
    public void setParentGroup(AtomGroup parent) {parentGroup = parent;}

    /**
     * Simulation in which this atom resides
     */
    public Simulation parentSimulation() {return parentSpecies().parentSimulation();}        
    /**
     * Phase in which this atom resides
     */
    public Phase parentPhase() {return parentGroup.parentPhase();}

    public Species parentSpecies() {return parentSpeciesAgent().parentSpecies();}
    
    public SpeciesAgent parentSpeciesAgent() {return parentGroup.parentSpeciesAgent();}
/*   linked list of bonds
    Bond firstBond;
    
    */
    public void sendToReservoir() {creator().reservoir().addAtom(this);}
    public AtomFactory creator() {return type.creator();}
    
    /**
     * Coordinates of this atom.
     * When the atom is constructed the coordinate class is provided by the 
     * governing Space for the simulation.
     */
    public final Space.Coordinate coord;
            
    /**
     * Integer assigned to this atom by its parent molecule.
     * Assigned during construction of atom.
     */
    public final int index() {return atomIndex;}
    public final void setIndex(int i) {atomIndex = i;}
        
    /**
     * Sets atom following this one in linked list, and sets this to be that atom's previous atom in list
     * 
     * @param atom the next atom in the list
     */
    public final void setNextAtom(Atom atom) {
        coord.setNextAtom(atom);
        /*
        nextAtom = atom;
        if(atom != null) atom.previousAtom = this;
        */
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
    
    public int depth() {return parentGroup.depth()+1;}
    
    public final boolean preceeds(Atom atom) {
        //want to return true if atoms are the same atoms
        if(atom == null) return true;
        if(this.parentGroup == atom.parentGroup()) return this.index() <= atom.index();//works also if both parentGroups are null
        int thisDepth = depth();
        int atomDepth = atom.depth();
        if(thisDepth == atomDepth) return this.parentGroup.preceeds(atom.parentGroup());
        else if(thisDepth < atomDepth) return this.preceeds(atom.parentGroup());
        else /*if(this.depth > atom.depth)*/ return this.parentGroup.preceeds(atom);
    }

    public final boolean isDescendedFrom(AtomGroup group) {
        return (this == group) || (parentGroup != null && parentGroup.isDescendedFrom(group));
    }
     /*   AtomGroup ancestor = parentGroup;
        while(ancestor != null) {
            if(ancestor == group) return true;
            ancestor = ancestor.parentGroup();
        }
        return false;
    }*/
    
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

    public Integrator.Agent ia;
        
    /**
    * Color of the atom when drawn on the screen
    * This color is set by the colorScheme object in the atom's species
    */
    Color color = Default.ATOM_COLOR;
        
    /**
    * Identifier of atom within molecule.
    * Assigned by parent molecule when invoking Atom constructor.
    */
     int atomIndex;
    
    private Atom nextAtom, previousAtom;
        
    public final AtomType type;
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