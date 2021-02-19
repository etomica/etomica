/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.dcvgcmd;

import etomica.action.AtomActionRandomizeVelocity;
import etomica.box.Box;
import etomica.integrator.IntegratorBoxFasterer;
import etomica.integrator.mcmove.MCMoveInsertDeleteFasterer;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeArrayList;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.util.random.IRandom;

/**
 * @author kofke
 * <p>
 * Extends insert/delete mcmove class by permitting insertion/deletion in a
 * subregion defined by a range of values of the z coordinate.
 */
public class MyMCMoveFasterer extends MCMoveInsertDeleteFasterer {

	public MyMCMoveFasterer(IntegratorBoxFasterer integrator, IRandom random,
							Space space, double zFraction) {
		super(integrator.getPotentialCompute(), random, space);
		position = space.makeVector();
		setZFraction(zFraction);
		this.integrator = integrator;
		randomizer = new AtomActionRandomizeVelocity(0, random);
		activeMolecules = new MoleculeArrayList();
	}

	public void setBox(Box p) {
		super.setBox(p);
	}

	/**
	 * Chooses and performs with equal probability an elementary molecule insertion
	 * or deletion.
	 */
	public boolean doTrial() {
		insert = (random.nextInt(2) == 0);
		if (insert) {
			uOld = 0.0;
			if (!reservoir.isEmpty()) {
				testMolecule = reservoir.remove(reservoir.size() - 1);
			} else {
				testMolecule = species.makeMolecule();
			}

			position.E(positionSource.randomPosition());
			double z = position.getX(2);
			double zBoundary = box.getBoundary().getBoxSize().getX(2);
			z += 0.5 * zBoundary;
			z *= zFraction;
			if (leftSide) {
				z = -0.5 * zBoundary + z;
			} else {
				z = 0.5 * zBoundary - z;
			}
			position.setX(2, z); //multiply z-coordinate by zFraction
			atomTranslator.setDestination(position);
			atomTranslator.actionPerformed(testMolecule);

			box.addMolecule(testMolecule);
		} else {//delete
			if (activeMolecules.size() == 0) {
				testMolecule = null;
				return false;
			}
			testMoleculeIndex = random.nextInt(activeMolecules.size());
			testMolecule = activeMolecules.get(testMoleculeIndex);
			uOld = potentialMaster.computeOneOldMolecule(testMolecule);
		}
		uNew = Double.NaN;
		return true;
	}//end of doTrial

	public double getChi(double temperature) {//note that moleculeCount() gives the number of molecules after the trial is attempted
		if (insert) {
			uNew = potentialMaster.computeOneMolecule(testMolecule);
		} else {
			uNew = 0;
		}
		double b = uOld - uNew;
		if (insert) b += mu;
		else b -= mu;

		double a = insert ? zFraction * box.getBoundary().volume() / (activeMolecules.size() + 1)
				: activeMolecules.size() / zFraction / box.getBoundary().volume();
		return a * Math.exp(b / temperature);
	}

	public void acceptNotify() {
		super.acceptNotify();
		if (!insert) {
			activeMolecules.remove(testMoleculeIndex);
			deltaN--;
		} else {
			activeMolecules.add(testMolecule);
			randomizer.setTemperature(integrator.getTemperature());
			randomizer.actionPerformed(testMolecule.getChildList().get(0));
			deltaN++;
		}
	}

	public void setupActiveAtoms() {
		activeMolecules.clear();
		double zBoundary = box.getBoundary().getBoxSize().getX(2);
		double zmin = leftSide ? -0.5 * zBoundary : (0.5 - zFraction) * zBoundary;
		double zmax = zmin + zFraction * zBoundary;
		int nMolecules = moleculeList.size();
		for (int i = 0; i < nMolecules; i++) {
			IMolecule molecule = moleculeList.get(i);

			double z = molecule.getChildList().get(0).getPosition().getX(2);
			if (z < zmin || z > zmax) continue;
			activeMolecules.add(molecule);
		}
	}

	public int getDeltaN() {
		return deltaN;
	}

	private double zFraction;
	private int deltaN = 0;
	private Vector position;
	private boolean leftSide;
	private final MoleculeArrayList activeMolecules;
	private IMoleculeList moleculeList;
	private final AtomActionRandomizeVelocity randomizer;
	private final IntegratorBoxFasterer integrator;
	protected int testMoleculeIndex;

	/**
	 * Returns the zFraction.
	 *
	 * @return double
	 */
	public double getZFraction() {
		return leftSide ? -zFraction : +zFraction;
	}

	/**
	 * Sets the zFraction, the fraction of the box volume into which atoms are
	 * inserted or deleted.  The volume is on the far left or right side (in
	 * the z dimension).  To put the volume on the left side (negative z),
	 * specify a negative z fraction (a positive value will result in the
	 * volume being on the positive z side.
	 */
	public void setZFraction(double zFraction) {
		this.zFraction = zFraction;
		leftSide = zFraction < 0.0;
		this.zFraction = Math.abs(zFraction);
	}

	public void setSpecies(ISpecies s) {
		super.setSpecies(s);
		moleculeList = box.getMoleculeList(s);
	}
}
