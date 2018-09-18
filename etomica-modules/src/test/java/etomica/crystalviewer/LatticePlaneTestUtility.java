/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.crystalviewer;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.lattice.*;
import etomica.lattice.crystal.*;
import etomica.simulation.Simulation;
import etomica.space.BoundaryDeformableLattice;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesSpheresMono;
import org.junit.jupiter.api.Assertions;

public class LatticePlaneTestUtility {

	public final int SIMPLE_CUBIC  = 0;
	public final int TETRAGONAL    = 1;
	public final int HEXAGONAL     = 2;
	public final int ORTHORHOMBIC  = 3;
	public final int MONOCLINIC    = 4;
	public final int TRICLINIC     = 5;
	public final int FCC           = 6;
	public final int BCC           = 7;
	public final int HCP           = 8;
	public final int CUBIC_DIAMOND = 9;

	private Simulation sim = null;
	private ISpecies species = null;
	private BravaisLattice lattice = null;
	private Box box = null;
	private LatticePlane latticePlane = null;
	private final Space3D space;

	public LatticePlaneTestUtility() {
    	// Create a simulation
		this.space = Space3D.getInstance();
    	sim = new Simulation(space);
    	space.makeVector();

    	// Create a species with one atom / molecule
    	species = new SpeciesSpheresMono(sim, space);
    	sim.addSpecies(species);

	}

	public void createLatticeAndBox(int index, int[] millerIndices, int[] boxSize) {

		BasisMonatomic basisMonatomic = null;

		// Create the lattice
		switch (index) {

			case SIMPLE_CUBIC:
				lattice = new LatticeCubicSimple(space);
				break;
			case TETRAGONAL:
				basisMonatomic = new BasisMonatomic(space);
				lattice = new BravaisLatticeCrystal(new PrimitiveTetragonal(space), basisMonatomic);
				break;
			case HEXAGONAL:
				basisMonatomic = new BasisMonatomic(space);
				lattice = new BravaisLatticeCrystal(new PrimitiveHexagonal(space), basisMonatomic);
				break;
			case ORTHORHOMBIC:
				basisMonatomic = new BasisMonatomic(space);
				lattice = new BravaisLatticeCrystal(new PrimitiveOrthorhombic(space), basisMonatomic);
				break;
			case MONOCLINIC:
				basisMonatomic = new BasisMonatomic(space);
				lattice = new BravaisLatticeCrystal(new PrimitiveMonoclinic(space), basisMonatomic);
				break;
			case TRICLINIC:
				basisMonatomic = new BasisMonatomic(space);
				lattice = new BravaisLatticeCrystal(new PrimitiveTriclinic(space), basisMonatomic);
				break;
			case FCC:
				lattice = new LatticeCubicFcc(space);
				break;
			case BCC:
				lattice = new LatticeCubicBcc(space);
				break;
			case HCP:
				lattice = new LatticeHcp(space);
				break;
			case CUBIC_DIAMOND:
				lattice = new LatticeCubicDiamond(space);
				break;
			default:
				break;
		}

		// Create the lattice plane
		latticePlane =
				new LatticePlane(lattice.getPrimitive(), millerIndices);
		setLatticePlanePosition(0.0);

		// Create a box
		if (box != null) {
			sim.removeBox(box);
		}
		box = sim.makeBox(new BoundaryDeformableLattice(lattice.getPrimitive(), boxSize));

	}

	public void setDimensions(int size) {
		((BoundaryDeformableLattice) box.getBoundary()).updateSize(lattice.getPrimitive(), new int[]{size, size, size});

		// Set the number of atoms
        int numAtoms = size*size*size;

	    // Set the number of molecules for the box and
	    // initialze the positions.
	    box.setNMolecules(species, numAtoms);
	    ConfigurationLattice config = new ConfigurationLattice(lattice, space);
	    config.initializeCoordinates(box);

	}

	public void setLatticePlanePosition(double pos) {
		latticePlane.setPosition(pos);
	}

	public double getLatticePlaneSpacePosition() {
		return latticePlane.getSpacePosition();
	}

	public Simulation getSimulation() {
		return sim;
	}

	public ISpecies getSpecies() {
		return species;
	}

	public BravaisLattice getLattice() {
		return lattice;
	}

	public Box getBox() {
		return box;
	}

	public LatticePlane getLatticePlane() {
		return latticePlane;
	}

	public void assertInPlane(IAtom a) {
		Assertions.assertTrue(this.getLatticePlane().inPlane(a.getPosition()), "Atom position " + a.getPosition() + " should be in plane but is not");
	}

	public void assertNotInPlane(IAtom a) {
		Assertions.assertFalse(this.getLatticePlane().inPlane(a.getPosition()), "Atom position " + a.getPosition() + " should not be in plane but is");
	}
}
