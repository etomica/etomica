/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode.nptdemo;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.graphics.ColorSchemeCollectiveAgent;
import etomica.normalmode.CoordinateDefinition;
import etomica.potential.compute.NeighborIterator;
import etomica.potential.compute.NeighborManager;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;

import java.awt.*;

/**
 * Color atoms based on being neighbors of the reference atom
 *
 * @author Andrew Schultz
 */
public class ColorSchemeScaledOverlap extends ColorSchemeCollectiveAgent {

    public ColorSchemeScaledOverlap(Space space, NeighborManager neighborManager, CoordinateDefinition coordinateDefinition) {
        super(coordinateDefinition.getBox());
        this.coordinateDefinition = coordinateDefinition;
        Box box = coordinateDefinition.getBox();
        nOverlaps = new int[box.getLeafList().size()];
        neighborIterator = neighborManager.makeNeighborIterator();
        pi = space.makeVector();
        pj = space.makeVector();
        dr = space.makeVector();
    }

    public void setPressure(double newPressure) {
        pressure = newPressure;
    }
    
    public double getPressure() {
        return pressure;
    }

    public void setDisplayDensity(double newDisplayDensity) {
        displayDensity = newDisplayDensity;
    }

    public double getDisplayDensity() {
        return displayDensity;
    }
    
    public synchronized void colorAllAtoms() {
        
        Box box = coordinateDefinition.getBox();
        IAtomList leafList = box.getLeafList();
        double vOld = box.getBoundary().volume();
        int nAtoms = box.getLeafList().size();
        double vNew = nAtoms/displayDensity;
        double rScale = Math.sqrt(vNew/vOld);
        double latticeScale = Math.exp((pressure*(vNew-vOld))/((nAtoms-1)*1*2))/rScale;
        // T=1, D=2
        
        
        double sigma = 1.0;
        double scaledSig = sigma/rScale;
        double sig2 = scaledSig*scaledSig;
        
        //color all atoms according to their type
        int nLeaf = leafList.size();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            nOverlaps[iLeaf] = 0;
        }
        Boundary boundary = box.getBoundary();
        for (int i = 0; i<leafList.size(); i++) {
            //color blue the neighbor atoms in same group
            IAtom atom = leafList.get(i);
            pi.E(atom.getPosition());
            Vector l = coordinateDefinition.getLatticePosition(atom);
            pi.ME(l);
            pi.TE(latticeScale);
            pi.PE(l);
            int finalI = i;

            neighborIterator.iterUpNeighbors(i, new NeighborIterator.NeighborConsumer() {
                @Override
                public void accept(IAtom jAtom, Vector rij, int n) {
                    pj.E(jAtom.getPosition());
                    Vector lj = coordinateDefinition.getLatticePosition(jAtom);
                    pj.ME(lj);
                    pj.TE(latticeScale);
                    pj.PE(lj);

                    dr.Ev1Mv2(pi, pj);
                    boundary.nearestImage(dr);
                    double r2 = dr.squared();
                    if (r2 < sig2) {
                        nOverlaps[finalI]++;
                        nOverlaps[jAtom.getLeafIndex()]++;
                    }
                }
            });
        }
        for (int i = 0; i<leafList.size(); i++) {
            //color green the target atom 
            agentManager.setAgent(leafList.get(i), colors[nOverlaps[i]]);
        }
    }

    private final NeighborIterator neighborIterator;
    protected final Vector dr;
    protected final Vector pi, pj;
    protected final int[] nOverlaps;
    protected double pressure, displayDensity;
    protected final CoordinateDefinition coordinateDefinition;
    protected Color[] colors = new Color[]{Color.RED, Color.BLUE, Color.GREEN, Color.BLACK, Color.CYAN, Color.PINK, new Color(0.5f, 0.0f, 0.5f)};
}
