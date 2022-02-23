package etomica.potential;

import etomica.atom.IAtom;

public interface BondingInfo {
    boolean skipBondedPair(boolean isPureAtoms, IAtom iAtom, IAtom jAtom);

    int n(boolean isPureAtoms, IAtom iAtom, IAtom jAtom);

    boolean isOnlyRigidMolecules();

    static BondingInfo noBonding() {
        return new BondingInfo() {
            @Override
            public boolean skipBondedPair(boolean isPureAtoms, IAtom iAtom, IAtom jAtom) {
                return !isPureAtoms && iAtom.getParentGroup() == jAtom.getParentGroup();
            }

            @Override
            public int n(boolean isPureAtoms, IAtom iAtom, IAtom jAtom) {
                return 0;
            }

            @Override
            public boolean isOnlyRigidMolecules() {
                return true;
            }
        };
    }
}
