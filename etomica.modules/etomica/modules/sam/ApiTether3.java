package etomica.modules.sam;

import etomica.action.AtomsetAction;
import etomica.api.IAtom;
import etomica.api.IAtomLeaf;
import etomica.api.IAtomSet;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.ISpecies;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomPair;
import etomica.atom.iterator.AtomsetIteratorPDT;
import etomica.atom.iterator.IteratorDirective.Direction;

/**
 * Returns first leaf atom of each polymer molecule and the atom its bonded to.
 */
public class ApiTether3 implements AtomsetIteratorPDT, AtomAgentManager.AgentSource {

    public ApiTether3(ISpecies polymerSpecies) {
        this.polymerSpecies = polymerSpecies;
        pair = new AtomPair();
    }
    
    public void setBox(IBox newBox) {
        if (box != newBox) {
            if (box != null) {
                bondManager.dispose();
            }
            bondManager = new AtomAgentManager(this, newBox);
        }
        box = newBox;
        polymerList = box.getMoleculeList(polymerSpecies);
    }

    public void setBondedSurfaceAtoms(IMolecule polymerMolecule, IAtomSet surfaceAtoms) {
        bondManager.setAgent(polymerMolecule, surfaceAtoms);
    }

    public void allAtoms(AtomsetAction action) {
        if (targetMolecule != null) {
            pair.atom0 = targetMolecule.getChildList().getAtom(0);
            for (int j=0; j<3; j++) {
                pair.atom1 = ((IAtomSet)(IAtom)bondManager.getAgent(targetMolecule)).getAtom(j);
                action.actionPerformed(pair);
            }
            return;
        }
        int nMolecules = polymerList.getAtomCount();
        for (int i=0; i<nMolecules; i++) {
            IMolecule molecule = (IMolecule)polymerList.getAtom(i);
            pair.atom0 = molecule.getChildList().getAtom(0);
            for (int j=0; j<3; j++) {
                pair.atom1 = ((IAtomSet)(IAtom)bondManager.getAgent(targetMolecule)).getAtom(j);
                action.actionPerformed(pair);
            }
        }
    }

    public int nBody() {
        return 2;
    }

    public IAtomSet next() {
        if (targetMolecule != null) {
            if (cursor == 0) {
                pair.atom0 = targetMolecule.getChildList().getAtom(0);
                pair.atom1 = ((IAtomSet)(IAtom)bondManager.getAgent(targetMolecule)).getAtom(surfaceCursor);
                surfaceCursor++;
                if (surfaceCursor == 3) {
                    cursor = 1;
                }
                return pair;
            }
            return null;
        }
        if (cursor > polymerList.getAtomCount() -1) {
            return null;
        }
        IMolecule molecule = (IMolecule)polymerList.getAtom(cursor);
        pair.atom0 = molecule.getChildList().getAtom(0);
        pair.atom1 = ((IAtomSet)bondManager.getAgent(molecule)).getAtom(surfaceCursor);
        surfaceCursor++;
        if (surfaceCursor == 3) {
            surfaceCursor = 0;
            cursor++;
        }
        return pair;
    }

    public void reset() {
        if (cursor > -1) {
            cursor = 0;
            surfaceCursor = 0;
        }
    }

    public int size() {
        return targetMolecule == null ? polymerList.getAtomCount()*3 : 3;
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

    public Class<IAtomSet> getAgentClass() {
        return IAtomSet.class;
    }

    public Object makeAgent(IAtom a) {
        return null;
    }

    public void releaseAgent(Object agent, IAtom atom) {
    }

    protected IBox box;
    protected final ISpecies polymerSpecies;
    protected IAtomSet polymerList;
    protected AtomAgentManager bondManager;
    protected int cursor;
    protected int surfaceCursor;
    protected IMolecule targetMolecule;
    protected final AtomPair pair;
}
