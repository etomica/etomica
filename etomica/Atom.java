package etomica;
import java.awt.Color;

 /**
  * Object corresponding to one physical atom or (in subclasses) group of atoms.
  * Each atom holds one instance of a Coordinate that is constructed by the governing
  * space class. 
  * all simulation kinetics and dynamics are performed by operating on the atom's coordinate.
  * In addition, an Atom has a type that holds information needed to draw the atom,
  * and possibly information to be used as parameters in a potential.  Each atom
  * holds an Integrator.Agent object that may be used to store information needed
  * by the integrator.
  * 
  * @author David Kofke
  * @author C. Daniel Barnes
  */
public class Atom implements java.io.Serializable {

    public static String getVersion() {return "Atom:01.08.08";}
    
    public Atom(Space space, AtomType t) {
        type = t;
        coord = space.makeCoordinate(this);//must follow setting of type field
//        coord.setMass(type.getMass());//handled in type.initialize statement
        if(atomLinkCount > 0) atomLink = new AtomLinker[atomLinkCount];//this is to be removed
        if(atomListCount > 0) atomList = new AtomList[atomListCount];
        for(int i=0; i<atomListCount; i++) {
            atomList[i] = new AtomList();
        }
        type.initialize(this);
    }
                    
    /**
     * Assigns the atom's integrator agent to the given instance.
     */
    public void setIntegratorAgent(Integrator.Agent ia) {this.ia = ia;}
            
    public final AtomGroup parentGroup() {return parentGroup;}
    
    public void setParentGroup(AtomGroup parent) {parentGroup = parent;}
    
    public int leafAtomCount() {return (type instanceof AtomType.Wall) ? 0 : 1;}
    
    public Atom firstChildAtom() {return this;}
    public Atom lastChildAtom() {return this;}
    public Atom firstLeafAtom() {return this;}
    public Atom lastLeafAtom() {return this;}

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

//   linked lists of bonds
    public BondLinker firstUpBond;
    public BondLinker firstDownBond;
    
    public void sendToReservoir() {
        if(parentGroup != null) parentGroup.removeAtom(this);
        creator().reservoir().addAtom(this);
    }
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
    public String signature() {return atomIndex + " " + parentGroup.signature();}
    public final String toString() {return "Atom(" + signature() + ")";}
        
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
    public final void clearPreviousAtom() {coord.clearPreviousAtom();}

    /**
     * @return the next atom in the linked list of atoms
     * @see setNextAtom
     */
    public final Atom nextAtom() {return coord.nextAtom();}

    /**
     * @return the previous atom in the linked list of atoms
     * @see setNextAtom
     */
    public final Atom previousAtom() {return coord.previousAtom();} 
    
    public int depth() {return (parentGroup != null) ? parentGroup.depth()+1 : 0;}
    
    public final boolean preceeds(Atom atom) {
        //want to return false if atoms are the same atoms
        if(atom == null) return true;
        if(this.parentGroup == atom.parentGroup()) return this.index() < atom.index();//works also if both parentGroups are null
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
    
//    private Atom nextAtom, previousAtom;
        
    public final AtomType type;
    protected AtomGroup parentGroup;
    
    public Object[] agents;
    
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
    public AtomList[] atomList;
    
    /**
     * Counter of the number of atom link index requests that have been fielded so far.
     * @deprecated
     */
    private static int atomLinkCount = 0;
    
    //replaces atomLinkCount
    private static int atomListCount = 0;
    
    /**
     * Returns a unique index that can be used by a simulation element to set up a linked
     * list of atoms associated with each atom.  If idx is the value returned by this
     * method, the element is permitted to set up a linked list of atoms beginning
     * with atomLink[idx].
     * @deprecated
     */
//    public static int requestAtomLinkIndex() {return atomLinkCount++;}
    
    //replaces AtomLinkIndex
    public static int requestAtomListIndex() {return atomListCount++;}
    
    
    /**
     * Interface for an object that makes an agent to be placed in each atom
     * upon construction.  AgentSource objects register with the AtomFactory
     * the produces the atom.
     */
    public interface AgentSource {
        public Object makeAgent(Atom a);
    }
    
}//end of Atom