package etomica.graphics;

import java.awt.Color;

import etomica.atom.Atom;
import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.nbr.PotentialMasterNbr;
import etomica.nbr.cell.Api1ACell;
import etomica.phase.Phase;
import etomica.simulation.Simulation;

public class ColorSchemeNeighbor extends ColorSchemeCollective {
    
    public ColorSchemeNeighbor(Simulation sim, Phase phase) {
        super(phase);
        typeColorScheme = new ColorSchemeByType();
        allIterator = new AtomIteratorLeafAtoms();
        allIterator.setPhase(phase);
        nbrIterator = new Api1ACell(sim.space().D(), sim.getDefaults().atomSize,((PotentialMasterNbr)sim.getPotentialMaster()).getCellAgentManager());
        nbrIterator.setDirection(null);
        nbrIterator.setPhase(phase);
    }
    
    public void colorAllAtoms() {
        Color[] atomColors = agentManager.getAgents();
        allIterator.reset();
		//color all atoms according to their type
        while(allIterator.hasNext()) {
            AtomLeaf atom = (AtomLeaf)allIterator.nextAtom();
            atomColors[atom.getGlobalIndex()] = typeColorScheme.getAtomColor(atom);//Color.green;
        }
        if (referenceAtom == null) {
            return;
        }
        //color blue the neighbor atoms in same group
        nbrIterator.reset();
        while(nbrIterator.hasNext()) {
            Atom  atom = nbrIterator.nextPair().atom1;
            if(atom.getType() == referenceAtom.getType()) {
                atomColors[atom.getGlobalIndex()] = Color.blue;
            } else {
                atomColors[atom.getGlobalIndex()] = Color.yellow;
            }
        }
        //color green the target atom 
        atomColors[referenceAtom.getGlobalIndex()] = Color.green;
    }
    
    public void setAtom(Atom a) {
        referenceAtom = a;
        nbrIterator.setTarget(a);
    }

    public Atom getAtom() {
        return referenceAtom;
    }
    
    private static final long serialVersionUID = 1L;
    private Atom referenceAtom;
    private final Api1ACell nbrIterator;
    private final AtomIteratorLeafAtoms allIterator;
    private final ColorSchemeByType typeColorScheme;
}