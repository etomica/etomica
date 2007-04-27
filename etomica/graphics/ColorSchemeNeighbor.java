package etomica.graphics;

import java.awt.Color;

import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomPair;
import etomica.atom.IAtom;
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
        Color[] atomColors = agentManager.getAgents();
		//color all atoms according to their type
        int nLeaf = leafList.size();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            AtomLeaf atom = (AtomLeaf)leafList.get(iLeaf);
            atomColors[atom.getGlobalIndex()] = typeColorScheme.getAtomColor(atom);//Color.green;
        }
        if (referenceAtom == null) {
            return;
        }
        //color blue the neighbor atoms in same group
        nbrIterator.reset();
        for (AtomPair pair = nbrIterator.nextPair(); pair != null;
             pair = nbrIterator.nextPair()) {
            IAtom atom = pair.atom1;
            if(atom.getType() == referenceAtom.getType()) {
                atomColors[atom.getGlobalIndex()] = Color.blue;
            } else {
                atomColors[atom.getGlobalIndex()] = Color.yellow;
            }
        }
        //color green the target atom 
        atomColors[referenceAtom.getGlobalIndex()] = Color.green;
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