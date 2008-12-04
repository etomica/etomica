package etomica.modules.sam;

import etomica.action.AtomsetAction;
import etomica.api.IAtom;
import etomica.api.IAtomLeaf;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.atom.AtomPair;
import etomica.atom.MoleculeAgentManager;
import etomica.atom.iterator.AtomsetIteratorPDT;
import etomica.atom.iterator.IteratorDirective.Direction;

/**
 * Returns first leaf atom of each polymer molecule and the atom its bonded to.
 */
public class ApiTether implements AtomsetIteratorPDT, MoleculeAgentManager.MoleculeAgentSource {

    public ApiTether(ISimulation sim, ISpecies polymerSpecies) {
        this.sim = sim;
        this.polymerSpecies = polymerSpecies;
        pair = new AtomPair();
    }
    
    public void setBox(IBox newBox) {
        if (box != newBox) {
            if (box != null) {
                bondManager.dispose();
            }
            bondManager = new MoleculeAgentManager(sim, newBox, this);
        }
        box = newBox;
        polymerList = box.getMoleculeList(polymerSpecies);
    }

    public void setBondedSurfaceAtom(IMolecule polymerMolecule, IAtomLeaf surfaceAtom) {
        bondManager.setAgent(polymerMolecule, surfaceAtom);
    }

    public void allAtoms(AtomsetAction action) {
        if (targetMolecule != null) {
            pair.atom0 = targetMolecule.getChildList().getAtom(0);
            pair.atom1 = (IAtom)bondManager.getAgent(targetMolecule);
            action.actionPerformed(pair);
            return;
        }
        int nMolecules = polymerList.getAtomCount();
        for (int i=0; i<nMolecules; i++) {
            IMolecule molecule = (IMolecule)polymerList.getAtom(i);
            pair.atom0 = molecule.getChildList().getAtom(0);
            pair.atom1 = (IAtom)bondManager.getAgent(molecule);
            action.actionPerformed(pair);
        }
    }

    public int nBody() {
        return 2;
    }

    public IAtomList next() {
        if (targetMolecule != null) {
            if (cursor == 0) {
                pair.atom0 = targetMolecule.getChildList().getAtom(0);
                pair.atom1 = (IAtom)bondManager.getAgent(targetMolecule);
                cursor = 1;
                return pair;
            }
            return null;
        }
        if (cursor > polymerList.getAtomCount() -1) {
            return null;
        }
        IMolecule molecule = (IMolecule)polymerList.getAtom(cursor);
        pair.atom0 = molecule.getChildList().getAtom(0);
        cursor++;
        pair.atom1 = (IAtom)bondManager.getAgent(molecule);
        return pair;
    }

    public void reset() {
        if (cursor > -1) cursor = 0;
    }

    public int size() {
        return targetMolecule == null ? polymerList.getAtomCount() : 1;
    }

    public void unset() {
        cursor = Integer.MAX_VALUE;
    }

    public void setDirection(Direction direction) {
    }

    public void setTarget(IAtom targetAtom) {
        cursor = 0;
        if (targetAtom == null) {
            targetMolecule = null;
        }
        else if (targetAtom instanceof IMolecule) {
            targetMolecule = (IMolecule)targetAtom;
            if (targetMolecule.getType() != polymerSpecies) {
                //hmm, we shouldn't hit this case
                System.out.println("you make me sad");
                targetMolecule = null;
                cursor = -1;
            }
        }
        else {
            targetMolecule = ((IAtomLeaf)targetAtom).getParentGroup();
            if (targetMolecule.getType() != polymerSpecies) {
                //hmm, we shouldn't hit this case
                System.out.println("you make me sad");
                targetMolecule = null;
                cursor = -1;
            }
        }
    }

    public Class getMoleculeAgentClass() {
        return IAtomLeaf.class;
    }

    public Object makeAgent(IMolecule a) {
        return null;
    }

    public void releaseAgent(Object agent, IMolecule atom) {
    }

    protected final ISimulation sim;
    protected IBox box;
    protected final ISpecies polymerSpecies;
    protected IAtomList polymerList;
    protected MoleculeAgentManager bondManager;
    protected int cursor;
    protected IMolecule targetMolecule;
    protected final AtomPair pair;
}
