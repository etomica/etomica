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
    
    public Atom(Space space, AtomType type, AtomTreeNodeGroup parent) {
        this(space, type, AtomTreeNodeGroup.FACTORY, parent);
    }
    public Atom(Space space, AtomType type, 
                    AtomTreeNode.Factory nodeFactory, AtomTreeNodeGroup parent) {
        this(space, type, nodeFactory, IteratorFactorySimple.INSTANCE.atomSequencerFactory(), parent);
    }
    public Atom(Space space, AtomType type, 
                    AtomTreeNode.Factory nodeFactory,
                    AtomSequencer.Factory seqFactory, AtomTreeNodeGroup parent) {
        this.type = type;//do this first
        seq = seqFactory.makeSequencer(this);
        node = nodeFactory.makeNode(this, parent);//must follow setting of sequencer to permit addition to childlist of parent
                                                  //setting parent here in constructor bypasses the call to the parent's addAtomNotify method (which would be called if parent were set in setParent method)
        coord = space.makeCoordinate(this);//must follow setting of type field
        if(parent != null) parent.addAtomNotify(this);
        
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
        creator().reservoir().addAtom(this);
    }
    
    public AtomFactory creator() {return type.creator();}
    
    /**
     * Coordinates of this atom.
     * When the atom is constructed the coordinate class is provided by the 
     * governing Space for the simulation.
     */
    public final Space.Coordinate coord;
            
    public String signature() {return node.parentGroup().signature() + " " + node.index();}
    public final String toString() {return "Atom(" + signature() + ")";}    

    public Integrator.Agent ia;
                
//    private Atom nextAtom, previousAtom;

    public final AtomTreeNode node;
        
    public final AtomType type;
    
    public final AtomSequencer seq;// = new AtomSequencerSimple(this);
    
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