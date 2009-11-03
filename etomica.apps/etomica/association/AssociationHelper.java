package etomica.association;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.atom.AtomArrayList;

public class AssociationHelper {
    
    protected final AssociationManager associationManager;
    
    public AssociationHelper(AssociationManager associationManager) {
        this.associationManager = associationManager;
    }

    public int populateList(AtomArrayList mySmerList, IAtom atom){
        mySmerList.clear();
        mySmerList.add(atom);
        IAtomList bondList = associationManager.getAssociatedAtoms(atom);
        if (bondList.getAtomCount() > 2){
            return 0;
        }
        if (bondList.getAtomCount() == 0){
            return 1;
        }
        IAtom thisAtom = bondList.getAtom(0);
        mySmerList.add(thisAtom);
        IAtomList bondList1 = associationManager.getAssociatedAtoms(thisAtom);
        if (bondList1.getAtomCount() > 2){
            return 0;
        }
        IAtom previousAtom = atom;
        while (bondList1.getAtomCount() > 1){
            IAtom nextAtom = bondList1.getAtom(0);
            if (nextAtom == previousAtom){
                nextAtom = bondList1.getAtom(1);
            } 
            if (nextAtom == atom){
                return 1;
            }
            mySmerList.add(nextAtom);
            bondList1 = associationManager.getAssociatedAtoms(nextAtom);
            if (bondList1.getAtomCount() > 2){
                return 0;
            }
            previousAtom = thisAtom;
            thisAtom = nextAtom;
        }
        if (bondList.getAtomCount()>1){
            thisAtom = bondList.getAtom(1);
            mySmerList.add(thisAtom);
            bondList1 = associationManager.getAssociatedAtoms(thisAtom);
            if (bondList1.getAtomCount() > 2){
                return 0;
            }
            previousAtom = atom;
            while (bondList1.getAtomCount() > 1){
                IAtom nextAtom = bondList1.getAtom(0);
                if (nextAtom == previousAtom){
                    nextAtom = bondList1.getAtom(1);
                } 
                mySmerList.add(nextAtom);
                bondList1 = associationManager.getAssociatedAtoms(nextAtom);
                if (bondList1.getAtomCount() > 2){
                    return 0;
                }
                previousAtom = thisAtom;
                thisAtom = nextAtom;
            }
        }
        return 1;
    }

}
