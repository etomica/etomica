package etomica.graphics;

import java.awt.Color;

import etomica.atom.AtomSet;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.nbr.cell.Api1ACell;
import etomica.nbr.list.PotentialMasterList;
import etomica.box.Box;
import etomica.simulation.ISimulation;

public class ColorSchemeNeighbor extends ColorSchemeCollective {
    
    public ColorSchemeNeighbor(ISimulation sim, PotentialMasterList potentialMaster, Box box) {
        super(box);
        typeColorScheme = new ColorSchemeByType();
        leafList = box.getLeafList();
        nbrIterator = new Api1ACell(sim.getSpace().D(), 1.0, potentialMaster.getCellAgentManager());
        nbrIterator.setDirection(null);
        nbrIterator.setBox(box);
    }
    
    public void colorAllAtoms() {
		//color all atoms according to their type
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomPositioned atom = (IAtomPositioned)leafList.getAtom(iLeaf);
            agentManager.setAgent(atom, typeColorScheme.getAtomColor(atom));
        }
        if (referenceAtom == null) {
            return;
        }
        //color blue the neighbor atoms in same group
        nbrIterator.reset();
        for (AtomSet pair = nbrIterator.next(); pair != null;
             pair = nbrIterator.next()) {
            IAtom atom = pair.getAtom(1);
            if(atom.getType() == referenceAtom.getType()) {
                agentManager.setAgent(atom, Color.blue);
            } else {
                agentManager.setAgent(atom, Color.yellow);
            }
        }
        //color green the target atom 
        agentManager.setAgent(referenceAtom, Color.green);
    }
    
    public void setAtom(IAtom a) {
        referenceAtom = a;
        nbrIterator.setTarget(a);
    }

    public IAtom getAtom() {
        return referenceAtom;
    }
    
    private static final long serialVersionUID = 1L;
    private IAtom referenceAtom;
    private final Api1ACell nbrIterator;
    private final AtomSet leafList;
    private final ColorSchemeByType typeColorScheme;
}