/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.api.IAtomList;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.box.Box;
import etomica.config.ConfigurationLatticeSimple;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.ISpace;
import etomica.species.SpeciesSpheresMono;

public class BasisBigCell extends Basis {

    private static final long serialVersionUID = 1L;

    public BasisBigCell(ISpace space, Basis subBasis, int[] nSubCells) {
        super(makeScaledCoordinates(space, subBasis, nSubCells));
    }

    protected static IVector[] makeScaledCoordinates(ISpace space, Basis subBasis, int[] nSubCells) {
        // make pretend sim, species and box so we can find the appropriate coordinates
        Simulation sim = new Simulation(space);
        ISpecies species = new SpeciesSpheresMono(sim, space);
        sim.addSpecies(species);
        // we might be used in the context of a deformable boundary (non-rectangular primitive)
        // but because we only care about scaled coordinates, the deformation doesn't
        // change what our result should be.  so just pretend that it's rectangular.
        
        
        IBoundary boundary = new BoundaryRectangularPeriodic(space);
        Primitive primitive = new PrimitiveCubic(space);
        
        IBox box = new Box(boundary, space);
        sim.addBox(box);
        IVector vector = space.makeVector(nSubCells);
        box.getBoundary().setBoxSize(vector);
        int numMolecules = subBasis.getScaledCoordinates().length;
        for (int i=0; i<nSubCells.length; i++) {
            numMolecules *= nSubCells[i];
        }
        box.setNMolecules(species, numMolecules);
        ConfigurationLatticeSimple configLattice = new ConfigurationLatticeSimple(new BravaisLatticeCrystal(primitive, subBasis), space);
        configLattice.initializeCoordinates(box);

        IVector boxSize = boundary.getBoxSize();
        
        // retrieve real coordinates and scale them
        IAtomList atomList = box.getLeafList();
        IVectorMutable[] pos = new IVectorMutable[atomList.getAtomCount()];
        for (int i=0; i<atomList.getAtomCount(); i++) {
            pos[i] = space.makeVector();
            pos[i].E(atomList.getAtom(i).getPosition());
            pos[i].DE(boxSize);
            // coordinates now range from -0.5 to +0.5, we want 0.0 to 1.0
            pos[i].PE(0.5);
            
        }
        return pos;
    }
}
