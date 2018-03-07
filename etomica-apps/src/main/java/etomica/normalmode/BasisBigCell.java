/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.config.ConfigurationLatticeSimple;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.species.SpeciesSpheresMono;

public class BasisBigCell extends Basis {

    private static final long serialVersionUID = 1L;

    public BasisBigCell(Space space, Basis subBasis, int[] nSubCells) {
        super(makeScaledCoordinates(space, subBasis, nSubCells));
    }

    protected static Vector[] makeScaledCoordinates(Space space, Basis subBasis, int[] nSubCells) {
        // make pretend sim, species and box so we can find the appropriate coordinates
        Simulation sim = new Simulation(space);
        ISpecies species = new SpeciesSpheresMono(sim, space);
        sim.addSpecies(species);
        // we might be used in the context of a deformable boundary (non-rectangular primitive)
        // but because we only care about scaled coordinates, the deformation doesn't
        // change what our result should be.  so just pretend that it's rectangular.
        
        
        Boundary boundary = new BoundaryRectangularPeriodic(space);
        Primitive primitive = new PrimitiveCubic(space);
        
        Box box = new Box(boundary, space);
        sim.addBox(box);
        Vector vector = space.makeVector(nSubCells);
        box.getBoundary().setBoxSize(vector);
        int numMolecules = subBasis.getScaledCoordinates().length;
        for (int i=0; i<nSubCells.length; i++) {
            numMolecules *= nSubCells[i];
        }
        box.setNMolecules(species, numMolecules);
        ConfigurationLatticeSimple configLattice = new ConfigurationLatticeSimple(new BravaisLatticeCrystal(primitive, subBasis), space);
        configLattice.initializeCoordinates(box);

        Vector boxSize = boundary.getBoxSize();
        
        // retrieve real coordinates and scale them
        IAtomList atomList = box.getLeafList();
        Vector[] pos = new Vector[atomList.size()];
        for (int i = 0; i<atomList.size(); i++) {
            pos[i] = space.makeVector();
            pos[i].E(atomList.get(i).getPosition());
            pos[i].DE(boxSize);
            // coordinates now range from -0.5 to +0.5, we want 0.0 to 1.0
            pos[i].PE(0.5);
            
        }
        return pos;
    }
}
