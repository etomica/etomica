/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.dcvgcmd;

import etomica.action.AtomActionRandomizeVelocity;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.integrator.IntegratorBoxFasterer;
import etomica.integrator.mcmove.MCMoveInsertDeleteFasterer;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeArrayList;
import etomica.potential.PotentialMasterBonding;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.RotationTensor3D;
import etomica.species.ISpecies;
import etomica.util.random.IRandom;

/**
 * @author kofke
 * <p>
 * Extends insert/delete mcmove class by permitting insertion/deletion in a
 * subregion defined by a range of values of the z coordinate.
 */
public class MyMCMoveFasterer extends MCMoveInsertDeleteFasterer {

	public MyMCMoveFasterer(IntegratorBoxFasterer integrator, PotentialMasterBonding pmBonding, IRandom random,
							Space space, double zFraction) {
		super(null, random, space);
		this.pmBonding = pmBonding;
		position = space.makeVector();
		setZFraction(zFraction);
		this.integrator = integrator;
		randomizer = new AtomActionRandomizeVelocity(0, random);
		activeAtoms = new MoleculeArrayList();
	}

	public void setSpecies(ISpecies s, double temperature) {
		super.setSpecies(s);
		this.temperature = temperature;
		resevoirMolecule = species.makeMolecule();
	}

	/**
	 * Chooses and performs with equal probability an elementary molecule insertion
	 * or deletion.
	 */
	public boolean doTrial() {
		insert = (random.nextInt(2) == 0);
		if (insert) {
			uOld = 0.0;

			IAtomList atoms = resevoirMolecule.getChildList();
			double u = pmBonding.computeOneMolecule(resevoirMolecule);
			Vector dr = box.getSpace().makeVector();
			double stepSize = 0.1;
			for (int step = 0; step < 10; step++) {
				nReservoirTrials++;
				Vector r = atoms.get(random.nextInt(atoms.size() - 1)).getPosition();

				dr.setRandomCube(random);
				dr.TE(stepSize);
				r.PE(dr);
				double uNew = pmBonding.computeOneOldMolecule(resevoirMolecule);
				if (uNew > u && random.nextDouble() > Math.exp(-(uNew - u) / temperature)) {
					r.ME(dr);
					continue;
				}
				nReservoirTrialsAccepted++;
				if (nReservoirTrialsAccepted % 100000 == 0) {
					System.out.println(nReservoirTrials + " trials, " + (nReservoirTrialsAccepted / (double) nReservoirTrials) + " accepted");
				}
				u = uNew;
			}

			testMolecule = species.makeMolecule();
			for (int i = 0; i < atoms.size(); i++) {
				testMolecule.getChildList().get(i).getPosition().E(atoms.get(i).getPosition());
			}
			double u1 = random.nextDouble();
			double u2 = 2 * Math.PI * random.nextDouble();
			double u3 = 2 * Math.PI * random.nextDouble();
			double s1 = Math.sqrt(u1);
			double s2 = Math.sqrt(1 - u1);
			RotationTensor3D tensor = new RotationTensor3D();
			tensor.setQuaternions(new double[]{s1 * Math.sin(u2), s1 * Math.cos(u2),
					s2 * Math.sin(u3), s2 * Math.cos(u3)});
			for (IAtom a : testMolecule.getChildList()) {
				tensor.transform(a.getPosition());
			}

			position.E(positionSource.randomPosition());
			double z = position.getX(2);
			double zBoundary = box.getBoundary().getBoxSize().getX(2);
			z += 0.5 * zBoundary;
			z *= zFraction;
			z += zPadding;
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
			if (activeAtoms.size() == 0) {
				testMolecule = null;//added this line 09/19/02
				return false;
			}
			testMoleculeIndex = random.nextInt(activeAtoms.size());
			testMolecule = activeAtoms.get(testMoleculeIndex);
			uOld = integrator.getPotentialCompute().computeOneOldMolecule(testMolecule);
		}
		uNew = Double.NaN;
		return true;
	}//end of doTrial

	public double getChi(double temperature) {//note that moleculeCount() gives the number of molecules after the trial is attempted
		if (insert) {
			uNew = integrator.getPotentialCompute().computeOneMolecule(testMolecule);
//			System.out.println("uNew "+uNew);
		} else {
			uNew = 0;
		}
		double b = uOld - uNew;
		if (insert) b += mu;
		else b -= mu;

		double a = insert ? zFraction * box.getBoundary().volume() / (activeAtoms.size() + 1)
				: activeAtoms.size() / zFraction / box.getBoundary().volume();
		return a * Math.exp(b / temperature);
	}

	public void rejectNotify() {
		if (insert) {
			// rejected insertion - remove from box and return to reservoir
			box.removeMolecule(testMolecule);
			// test molecule is no longer in the simulation and should not be
			// returned by affectedAtoms
			testMolecule = null;
		}
	}

	public void acceptNotify() {
		if (insert) {
			integrator.getPotentialCompute().processAtomU(1);
			activeAtoms.add(testMolecule);
			randomizer.setTemperature(integrator.getTemperature());
			for (IAtom a : testMolecule.getChildList()) {
				randomizer.actionPerformed(a);
			}
			deltaN++;
		} else {
			integrator.getPotentialCompute().computeOneMolecule(testMolecule);
			integrator.getPotentialCompute().processAtomU(-1);
			// accepted deletion - remove from box and add to reservoir
			box.removeMolecule(testMolecule);
			// accepted deletion - remove from box and add to reservoir

//			System.out.println("removed "+testMolecule.getType().getIndex()+" "+testMoleculeIndex);
			activeAtoms.remove(testMoleculeIndex);
			deltaN--;
		}
	}

	public void setupActiveAtoms() {
		activeAtoms.clear();
		double zBoundary = box.getBoundary().getBoxSize().getX(2);
		double zmin = leftSide ? (zPadding + -0.5 * zBoundary) : (-zPadding + (0.5 - zFraction) * zBoundary);
		double zmax = zmin + zFraction * zBoundary;
		int nMolecules = moleculeList.size();
		for (int i = 0; i < nMolecules; i++) {
			IMolecule molecule = moleculeList.get(i);

			double z = molecule.getChildList().get(0).getPosition().getX(2);
			if (z < zmin || z > zmax) continue;
			activeAtoms.add(molecule);
		}
	}

	public int getDeltaN() {
		return deltaN;
	}

	private final PotentialMasterBonding pmBonding;
	private double zFraction, zPadding = 5;
	private int deltaN = 0;
	private final Vector position;
	private boolean leftSide;
	private final MoleculeArrayList activeAtoms;
	private IMoleculeList moleculeList;
	private final AtomActionRandomizeVelocity randomizer;
	private final IntegratorBoxFasterer integrator;
	protected int testMoleculeIndex;
	protected double temperature;
	protected IMolecule resevoirMolecule;
	public long nReservoirTrials, nReservoirTrialsAccepted;

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
