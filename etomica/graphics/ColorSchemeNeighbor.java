package etomica.graphics;

import java.awt.Color;

import etomica.api.IAtomLeaf;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.ISimulation;
import etomica.nbr.cell.Api1ACell;
import etomica.nbr.list.PotentialMasterList;

public class ColorSchemeNeighbor extends ColorSchemeCollectiveAgent {
    
    public ColorSchemeNeighbor(ISimulation sim, PotentialMasterList potentialMaster, IBox box, int dim) {
        super(box);
        typeColorScheme = new ColorSchemeByType(sim);
        leafList = box.getLeafList();
        nbrIterator = new Api1ACell(dim, 1.0, potentialMaster.getCellAgentManager());
        nbrIterator.setDirection(null);
        nbrIterator.setBox(box);
    }
    
    public void colorAllAtoms() {
		//color all atoms according to their type
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomLeaf atom = leafList.getAtom(iLeaf);
            agentManager.setAgent(atom, typeColorScheme.getAtomColor(atom));
        }
        if (referenceAtom == null) {
            return;
        }
        //color blue the neighbor atoms in same group
        nbrIterator.reset();
        for (IAtomList pair = nbrIterator.next(); pair != null;
             pair = nbrIterator.next()) {
            IAtomLeaf atom = pair.getAtom(1);
            if(atom.getType() == referenceAtom.getType()) {
                agentManager.setAgent(atom, Color.blue);
            } else {
                agentManager.setAgent(atom, Color.yellow);
            }
        }
        //color green the target atom 
        agentManager.setAgent(referenceAtom, Color.green);
    }
    
    public void setAtom(IAtomLeaf a) {
        referenceAtom = a;
        nbrIterator.setTarget(a);
    }

    public IAtomLeaf getAtom() {
        return referenceAtom;
    }
    
    private static final long serialVersionUID = 1L;
    private IAtomLeaf referenceAtom;
    private final Api1ACell nbrIterator;
    private final IAtomList leafList;
    private final ColorSchemeByType typeColorScheme;
}