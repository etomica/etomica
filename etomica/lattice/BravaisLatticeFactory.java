package etomica.lattice;

import etomica.*;

public class BravaisLatticeFactory extends AtomFactoryTree {
    
    public Atom build() {
        Atom group = new BravaisLattice(space, groupType);
        return build(group);
    }
}