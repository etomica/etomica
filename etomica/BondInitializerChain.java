package etomica;

/**
 * Adds bonds to an atom group such that atoms adjacent in the sequence
 * are bonded, thus forming the atoms into a chain.
 *
 * @author David Kofke
 */

public class BondInitializerChain extends BondInitializer {
    
    private AtomIteratorSequential iterator = new AtomIteratorSequential();
    private Bond bondType;
    
    public BondInitializerChain() {
        this(Bond.instance);
    }
    public BondInitializerChain(Bond type) {
        bondType = type;
    }
    
    public void makeBonds(Atom a) {
        if(!(a instanceof AtomGroup)) return;
        Atom a1 = null;
        Atom a2 = null;
        iterator.setBasis(a);
        iterator.reset();
        while(iterator.hasNext()) {
            a2 = iterator.next();
            if(a1 != null && a2 != null) bondType.makeBond(a1, a2);
            a1 = a2;
        }
    }
    
}//end of BondInitializerChain
