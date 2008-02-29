package etomica.atom;

import etomica.api.IBox;
import etomica.chem.elements.ElementSimple;

 /**
  * Object corresponding to one physical atom or group of atoms. Each atom holds
  * the following publicly accessible fields:
  * <ul>
  * <li>an AtomType instance (fieldname: type) that holds information this atom
  * has in common with other atoms made by the same factory
  * </ul>
  * <p>
  * @author David Kofke, Andrew Schultz, and C. Daniel Barnes
  * 
  */
public abstract class Atom implements IAtom, java.io.Serializable {

    public Atom(AtomType type) {
        this.type = type;
    }
    
    /**
     * Makes a simple atom.  Node is for a leaf atom; 
     * type is a sphere with unit mass and unit size, unique to the new atom; 
     * depth is 0.
     */
    public Atom() {
        this(makeAtomTypeSphere());                        
        setIndex(++INSTANCE_COUNT);//default index; changed when added to parent after construction
    }
    
    /**
     * Method to return a dummy AtomType that's valid (has parent explicitly set to null)
     */
    private static AtomTypeSphere makeAtomTypeSphere() {
        return new AtomTypeSphere(new ElementSimple("Simple",1),1);
    }
    
    public abstract String signature();
    
    public void setGlobalIndex(IBox box) {
        globalIndex = box.requestGlobalIndex();
    }

    public final int getGlobalIndex() {
        return globalIndex;
    }
    
    /**
     * @return the Atom type, holding properties held in common with other 
     * atoms made by this atom's factory.
     */
    public final AtomType getType() {
        return type;
    }

    public final void setIndex(int newIndex) {
        index = newIndex;
    }
    
    public final int getIndex() {
        return index;
    }
    
    protected final AtomType type;
    
    /**
     * Counter for number of times an atom is instantiated without a parent.  Used
     * to assign a unique index to such atoms.
     */
    private static int INSTANCE_COUNT = 0;
    
    private int globalIndex = -1;
    protected int index;
}
