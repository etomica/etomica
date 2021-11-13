/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.compute.NeighborIterator;
import etomica.potential.compute.NeighborManager;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;

import java.awt.*;

/**
 * Color atoms based on how many overlaps they have with their neighbors
 *
 * @author Andrew Schultz
 */
public class ColorSchemeOverlap extends ColorSchemeCollectiveAgent {

    public ColorSchemeOverlap(Space space, NeighborManager neighborManager, Box box) {
        super(box);
        leafList = box.getLeafList();
        nOverlaps = new int[leafList.size()];
        this.neighborIterator = neighborManager.makeNeighborIterator();
        dr = space.makeVector();
        boundary = box.getBoundary();
    }

    public synchronized void colorAllAtoms() {
        int nLeaf = leafList.size();
        if (nOverlaps.length != nLeaf) {
            nOverlaps = new int[nLeaf];
        }
        else {
            for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
                nOverlaps[iLeaf] = 0;
            }
        }
        for (int i = 0; i<leafList.size(); i++) {
            IAtom atom = leafList.get(i);
            int finalI = i;
            neighborIterator.iterUpNeighbors(i, new NeighborIterator.NeighborConsumer() {
                @Override
                public void accept(IAtom jAtom, Vector rij) {
                    double r2 = rij.squared();
                    if (r2 < sig2) {
                        // count overlaps for both i and j
                        nOverlaps[finalI]++;
                        nOverlaps[jAtom.getLeafIndex()]++;
                    }
                }
            });
        }
        for (int i = 0; i<leafList.size(); i++) {
            // set appropriate color for the # of overlaps for each atom
            agentManager.setAgent(leafList.get(i), colors[nOverlaps[i]]);
        }
    }

    public void setSigma(double newSigma) {
        sig2 = newSigma*newSigma;
    }

    protected final NeighborIterator neighborIterator;
    protected final Boundary boundary;
    protected final IAtomList leafList;
    protected final Vector dr;
    protected double sig2;
    protected int[] nOverlaps;
    protected Color[] colors = new Color[]{Color.RED, Color.BLUE, Color.GREEN, Color.BLACK, Color.CYAN, Color.PINK, new Color(0.5f, 0.0f, 0.5f)};
}
