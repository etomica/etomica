package etomica;

/**
 * Adds bonds to an atom group such that atoms adjacent in the sequence
 * are bonded, thus forming the atoms into a chain.
 *
 * @author David Kofke
 */

public class BondInitializerChain extends BondInitializer {
    
    private final AtomIteratorList iterator = new AtomIteratorList();
    private Bond bondType;
    
    public BondInitializerChain() {
        this(Bond.instance);
    }
    public BondInitializerChain(Bond type) {
        bondType = type;
    }
    
    public void makeBonds(Atom a) {
        if(a.node.isLeaf()) return;
        Atom a1 = null;
        Atom a2 = null;
        iterator.setList(((AtomTreeNodeGroup)a.node).childList);
        iterator.reset();
        while(iterator.hasNext()) {
            a2 = iterator.nextAtom();
            if(a1 != null && a2 != null) Bond.makeBond(a1, a2);
            a1 = a2;
        }
    }
    
}//end of BondInitializerChain
