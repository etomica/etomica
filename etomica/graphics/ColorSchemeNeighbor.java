package etomica.graphics;

import java.awt.Color;

import etomica.atom.Atom;
import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.nbr.PotentialMasterNbr;
import etomica.nbr.cell.Api1ACell;
import etomica.phase.Phase;
import etomica.simulation.Simulation;



/*
 * History
 * Created on Aug 23, 2005 by kofke
 */
public class ColorSchemeNeighbor extends ColorSchemeCollective {
    
    private Atom referenceAtom;
    private final Api1ACell nbrIterator;
    private final AtomIteratorLeafAtoms allIterator = new AtomIteratorLeafAtoms();
    private final ColorSchemeByType typeColorScheme = new ColorSchemeByType();
    
    public ColorSchemeNeighbor(Simulation sim) {
        super(sim);
        nbrIterator = new Api1ACell(sim.space().D(), sim.getDefaults().atomSize,((PotentialMasterNbr)sim.potentialMaster).getCellAgentManager());
        nbrIterator.setDirection(null);
    }
    
    public void colorAllAtoms(Phase phase) {
        allIterator.setPhase(phase);
        allIterator.reset();
		//color all atoms according to their type
        while(allIterator.hasNext()) {
            AtomLeaf atom = (AtomLeaf)allIterator.nextAtom();
            atomColors[atom.getGlobalIndex()] = typeColorScheme.getAtomColor(atom);//Color.green;
        }
        //color blue the neighbor atoms in same group
        nbrIterator.setPhase(phase);
        nbrIterator.reset();
        while(nbrIterator.hasNext()) {
            Atom  atom = nbrIterator.nextPair().atom1;
            if(atom.type == referenceAtom.type) {
                atomColors[atom.getGlobalIndex()] = Color.blue;
            } else {
                atomColors[atom.getGlobalIndex()] = Color.yellow;
            }
        }
        //color red the target atom 
        atomColors[referenceAtom.getGlobalIndex()] = Color.red;
    }
    
    public void setAtom(Atom a) {
        referenceAtom = a;
        nbrIterator.setTarget(a);
    }

    public Atom getAtom() {return referenceAtom;}
}