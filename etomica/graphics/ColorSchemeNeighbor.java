package etomica.graphics;

import java.awt.Color;

import etomica.atom.AtomArrayList;
import etomica.atom.AtomSet;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.nbr.PotentialMasterNbr;
import etomica.nbr.cell.Api1ACell;
import etomica.phase.Phase;
import etomica.simulation.Simulation;

public class ColorSchemeNeighbor extends ColorSchemeCollective {
    
    public ColorSchemeNeighbor(Simulation sim, Phase phase) {
        super(phase);
        typeColorScheme = new ColorSchemeByType();
        leafList = phase.getSpeciesMaster().getLeafList();
        nbrIterator = new Api1ACell(sim.getSpace().D(), sim.getDefaults().atomSize,((PotentialMasterNbr)sim.getPotentialMaster()).getCellAgentManager());
        nbrIterator.setDirection(null);
        nbrIterator.setPhase(phase);
    }
    
    public void colorAllAtoms() {
		//color all atoms according to their type
        int nLeaf = leafList.size();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomPositioned atom = (IAtomPositioned)leafList.get(iLeaf);
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
    private final AtomArrayList leafList;
    private final ColorSchemeByType typeColorScheme;
}