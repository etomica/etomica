package etomica.atom;

import java.io.Serializable;

import etomica.api.IMolecule;
import etomica.api.IMoleculeList;

/**
 * @author Tai Boon Tan  
 * 
 * 
 */
public class MoleculeToMoleculeSetFixed implements MoleculeToMoleculeList, MoleculeToIndex, Serializable {

    private static final long serialVersionUID = 1L;

    public MoleculeToMoleculeSetFixed() {
        moleculeArrayList = new MoleculeArrayList();
    }
    
    public void setArrayList(MoleculeArrayList list) {
        moleculeArrayList = list;
    }
    
    public IMoleculeList getMoleculeList(IMolecule molecule) {
        return moleculeArrayList;
    }
    
    public int getIndex(IMolecule molecule) {
        return moleculeArrayList.indexOf(molecule);
    }

    private MoleculeArrayList moleculeArrayList;
}
