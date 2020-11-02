package etomica.potential;

import etomica.atom.IAtom;

public interface BondingInfo {
    boolean skipBondedPair(boolean isPureAtoms, IAtom iAtom, IAtom jAtom);

    boolean isOnlyRigidMolecules();

    static BondingInfo noBonding() {
        return new BondingInfo() {
            @Override
            public boolean skipBondedPair(boolean isPureAtoms, IAtom iAtom, IAtom jAtom) {
                return !isPureAtoms && iAtom.getParentGroup() == jAtom.getParentGroup();
            }

            @Override
            public boolean isOnlyRigidMolecules() {
                return true;
            }
        };
    }

    static BondingInfo makeBondingInfo(PotentialMasterBonding pm) {
        if (pm != null) {
            return new BondingInfo() {
                @Override
                public boolean skipBondedPair(boolean isPureAtoms, IAtom iAtom, IAtom jAtom) {
                    if (!isPureAtoms && iAtom.getParentGroup() == jAtom.getParentGroup()) {
                        // ensure i < j
                        if (pm.isOnlyRigidMolecules) {
                            return true;
                        }
                        if (iAtom.getLeafIndex() > jAtom.getLeafIndex()) {
                            IAtom tmp = iAtom;
                            iAtom = jAtom;
                            jAtom = tmp;
                        }
                        int species = iAtom.getParentGroup().getType().getIndex();
                        int iChildIndex = iAtom.getIndex();
                        int jChildIndex = jAtom.getIndex();
                        int[] iBondedAtoms = pm.bondedAtoms[species][iChildIndex];
                        for (int iBondedAtom : iBondedAtoms) {
                            if (iBondedAtom == jChildIndex) return true;
                        }
                    }
                    return false;
                }

                @Override
                public boolean isOnlyRigidMolecules() {
                    return pm.isOnlyRigidMolecules;
                }
            };
        } else {
            return noBonding();
        }
    }
}
