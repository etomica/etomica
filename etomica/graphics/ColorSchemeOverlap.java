package etomica.graphics;

import java.awt.Color;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.space.ISpace;

/**
 * Color atoms based on how many overlaps they have with their neighbors
 *
 * @author Andrew Schultz
 */
public class ColorSchemeOverlap extends ColorSchemeCollectiveAgent {
    
    public ColorSchemeOverlap(ISpace space, PotentialMasterList potentialMaster, IBox box) {
        super(box);
        leafList = box.getLeafList();
        nOverlaps = new int[leafList.getAtomCount()];
        neighborManager = potentialMaster.getNeighborManager(box);
        dr = space.makeVector();
    }

    public synchronized void colorAllAtoms() {
        int nLeaf = leafList.getAtomCount();
        if (nOverlaps.length != nLeaf) {
            nOverlaps = new int[nLeaf];
        }
        else {
            for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
                nOverlaps[iLeaf] = 0;
            }
        }
        for (int i=0; i<leafList.getAtomCount(); i++) {
            IAtom atom = leafList.getAtom(i);
            IAtomList list = neighborManager.getDownList(atom)[0];
            IVector p = atom.getPosition();
            for (int j=0; j<list.getAtomCount(); j++) {
                IAtom jAtom = list.getAtom(j);
                double r2 = p.Mv1Squared(jAtom.getPosition());
                if (r2 < sig2) {
                    // count overlaps for both i and j
                    nOverlaps[i]++;
                    nOverlaps[jAtom.getLeafIndex()]++;
                }
            }
        }
        for (int i=0; i<leafList.getAtomCount(); i++) {
            // set appropriate color for the # of overlaps for each atom
            agentManager.setAgent(leafList.getAtom(i), colors[nOverlaps[i]]);
        }
    }

    public void setSigma(double newSigma) {
        sig2 = newSigma*newSigma;
    }

    private static final long serialVersionUID = 1L;
    private final NeighborListManager neighborManager;
    private final IAtomList leafList;
    protected final IVectorMutable dr;
    protected double sig2;
    protected int[] nOverlaps;
    protected Color[] colors = new Color[]{Color.RED, Color.BLUE, Color.GREEN, Color.BLACK, Color.CYAN, Color.PINK, new Color(0.5f, 0.0f, 0.5f)};
}