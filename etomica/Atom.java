package etomica;

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
        this(space, t, new AtomTreeNodeGroup());
    }
    public Atom(Space space, AtomType t, AtomTreeNode node) {
        type = t;
        coord = space.makeCoordinate(this);//must follow setting of type field
        this.node = node;
        node.setAtom(this);
        
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
            

//   linked lists of bonds
    public BondLinker firstUpBond;
    public BondLinker firstDownBond;
    
    public void sendToReservoir() {
        if(node.parentGroup() != null) node.parentGroup().node.removeAtom(this);
        creator().reservoir().addAtom(this);
    }
    public AtomFactory creator() {return type.creator();}
    
    /**
     * Coordinates of this atom.
     * When the atom is constructed the coordinate class is provided by the 
     * governing Space for the simulation.
     */
    public final Space.Coordinate coord;
            
    public String signature() {return node.index() + " " + node.parentGroup().signature();}
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
    
    /**
     * Returns true if this atom preceeds the given atom in the atom sequence.
     * Returns false if the given atom is this atom, or (of course) if the
     * given atom instead preceeds this one.
     */
    public final boolean preceeds(Atom atom) {
        //want to return false if atoms are the same atoms
        if(atom == null) return true;
        if(this.node.parentGroup() == atom.node.parentGroup()) return this.node.index() < atom.node.index();//works also if both parentGroups are null
        int thisDepth = node.depth();
        int atomDepth = atom.node.depth();
        if(thisDepth == atomDepth) return this.node.parentGroup().preceeds(atom.node.parentGroup());
        else if(thisDepth < atomDepth) return this.preceeds(atom.node.parentGroup());
        else /*if(this.depth > atom.depth)*/ return this.node.parentGroup().preceeds(atom);
    }

    public Integrator.Agent ia;
                
//    private Atom nextAtom, previousAtom;

    public final AtomTreeNode node;
        
    public final AtomType type;
    
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
    
    //replaces AtomLinkIndex
    public static int requestAtomListIndex() {return atomListCount++;}
    
    /**
     * Returns a unique index that can be used by a simulation element to set up a linked
     * list of atoms associated with each atom.  If idx is the value returned by this
     * method, the element is permitted to set up a linked list of atoms beginning
     * with atomLink[idx].
     * @deprecated
     */
//    public static int requestAtomLinkIndex() {return atomLinkCount++;}
    
    
    /**
     * Interface for an object that makes an agent to be placed in each atom
     * upon construction.  AgentSource objects register with the AtomFactory
     * the produces the atom.
     */
    public interface AgentSource {
        public Object makeAgent(Atom a);
    }
    
}//end of Atom