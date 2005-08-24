package etomica.graphics;

import java.awt.Color;

import etomica.Phase;
import etomica.Simulation;
import etomica.atom.Atom;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.nbr.cell.Api1ACell;



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
        nbrIterator = new Api1ACell(sim.space().D());
        nbrIterator.setDirection(null);
    }
    
    public void colorAllAtoms(Phase phase) {
        allIterator.setPhase(phase);
        allIterator.reset();
		//color all atoms according to their type
        while(allIterator.hasNext()) {
            Atom atom = allIterator.nextAtom();
            atom.allatomAgents[agentIndex] = typeColorScheme.atomColor(atom);//Color.green;
        }
        //color blue the neighbor atoms in same group
        nbrIterator.setPhase(phase);
        nbrIterator.reset();
        while(nbrIterator.hasNext()) {
            Atom  atom = nbrIterator.nextPair().atom1;
            if(atom.type == referenceAtom.type) {
                atom.allatomAgents[agentIndex] = Color.blue;
            } else {
                atom.allatomAgents[agentIndex] = Color.yellow;
            }
        }
        //color red the target atom 
        referenceAtom.allatomAgents[agentIndex] = Color.red;
    }
    
    public void setAtom(Atom a) {
        referenceAtom = a;
        nbrIterator.setTarget(a);
    }

    public Atom getAtom() {return referenceAtom;}
}