/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.chainequilibrium;

import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.IRandom;
import etomica.api.IVectorMutable;
import etomica.config.ConfigurationLattice;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.SpaceLattice;
import etomica.space.ISpace;

/**
 * Configuration that puts atoms randomly on a lattice.
 * 
 * @author Andrew Schultz
 */
public class ConfigurationLatticeRandom extends ConfigurationLattice {

    public ConfigurationLatticeRandom(SpaceLattice lattice, ISpace space, IRandom random) {
        super(lattice, space);
        this.random = random;
    }

    public void initializeCoordinates(IBox box) {
        IMoleculeList moleculeList = box.getMoleculeList();
        int sumOfMolecules = moleculeList.getMoleculeCount();
        if (sumOfMolecules == 0) {
            return;
        }
        int basisSize = 1;
        if (lattice instanceof BravaisLatticeCrystal) {
            basisSize = ((BravaisLatticeCrystal)lattice).getBasis().getScaledCoordinates().length;
        }
        int nCells = (int) Math.ceil((double) sumOfMolecules
                / (double) basisSize);

        // determine scaled shape of simulation volume
        IVectorMutable shape = space.makeVector();
        shape.E(box.getBoundary().getBoxSize());
        shape.PE(-boundaryPadding);
        IVectorMutable latticeConstantV = space.makeVector(lattice.getLatticeConstants());
        shape.DE(latticeConstantV);

        // determine number of cells in each direction
        int[] latticeDimensions = calculateLatticeDimensions(nCells, shape);
        int nSites = basisSize;
        for (int i=0; i<latticeDimensions.length; i++) {
            nSites *= latticeDimensions[i];
        }
        if (indexIterator.getD() > latticeDimensions.length) {
            int[] iteratorDimensions = new int[latticeDimensions.length+1];
            System.arraycopy(latticeDimensions, 0, iteratorDimensions, 0,
                    latticeDimensions.length);
            iteratorDimensions[latticeDimensions.length] = basisSize;
            indexIterator.setSize(iteratorDimensions);
        }
        else {
            indexIterator.setSize(latticeDimensions);
        }

        // determine lattice constant
        IVectorMutable latticeScaling = space.makeVector();
        if (rescalingToFitVolume) {
            // in favorable situations, this should be approximately equal
            // to 1.0
            latticeScaling.E(box.getBoundary().getBoxSize());
            latticeScaling.PE(-boundaryPadding);
            latticeScaling.DE(latticeConstantV);
            latticeScaling.DE(space.makeVector(latticeDimensions));
        } else {
            latticeScaling.E(1.0);
        }

        // determine amount to shift lattice so it is centered in volume
        IVectorMutable offset = space.makeVector();
        offset.E(box.getBoundary().getBoxSize());
        IVectorMutable vectorOfMax = space.makeVector();
        IVectorMutable vectorOfMin = space.makeVector();
        IVectorMutable site = space.makeVector();
        vectorOfMax.E(Double.NEGATIVE_INFINITY);
        vectorOfMin.E(Double.POSITIVE_INFINITY);

        indexIterator.reset();

        while (indexIterator.hasNext()) {
            site.E((IVectorMutable) lattice.site(indexIterator.next()));
            site.TE(latticeScaling);
            for (int i=0; i<site.getD(); i++) {
                vectorOfMax.setX(i, Math.max(site.getX(i),vectorOfMax.getX(i)));
                vectorOfMin.setX(i, Math.min(site.getX(i),vectorOfMin.getX(i)));
            }
        }
        offset.Ev1Mv2(vectorOfMax, vectorOfMin);
        offset.TE(-0.5);
        offset.ME(vectorOfMin);

        myLat = new MyLattice(lattice, latticeScaling, offset);

        // Place molecules
        indexIterator.reset();
        double voidFrac = (nSites - sumOfMolecules)/((double)nSites);
        double voidSum = 0;
        int siteCount = 0;
        boolean[] done = new boolean[sumOfMolecules];
        for (int j=0; j<sumOfMolecules; j++) {
            int i;
            do {
                i = random.nextInt(sumOfMolecules);
            }
            while (done[i]);
            done[i] = true;
            IMolecule a = moleculeList.getMolecule(i);
            int[] ii = indexIterator.next();
            siteCount++;
            // add voidFrac for each /site/ (not molecule)
            voidSum += voidFrac;
            while (voidSum > 1.0) {
                // we've gone through enough sites that we should insert a void
                // now.  Subtract one, but still add voidFrac since we're still
                // advancing one site.
                voidSum += voidFrac - 1;
                ii = indexIterator.next();
                siteCount++;
            }
            // initialize coordinates of child atoms
            a.getType().initializeConformation(a);

            atomActionTranslateTo.setDestination((IVectorMutable)myLat.site(ii));
            atomActionTranslateTo.actionPerformed(a);
        }
        if (nSites - siteCount > Math.ceil(1.0/(1.0-voidFrac))) {
            // nSites - siteCount = 0 is ideal.
            // indexIterator.next() would throw if nSites < siteCount
            // nSites - siteCount = 1 will be typical for cases where the void distribution can't be perfect
            // so we just need to check for nSites - siteCount > 1
            // for very low occupancy lattices, we'll do worse.
            throw new RuntimeException("Failed to properly iterate through the lattice sites "+nSites+" "+siteCount);
        }
    }

    protected final IRandom random;
}
