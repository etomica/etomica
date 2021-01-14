/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/*
 * Created on Apr 12, 2004
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package etomica.dcvgcmd;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorMD;
import etomica.integrator.mcmove.MCMoveManager;
import etomica.molecule.IMolecule;
import etomica.nbr.PotentialMasterHybrid;
import etomica.potential.PotentialMaster;
import etomica.species.ISpecies;
import etomica.units.Kelvin;
import etomica.util.random.IRandom;


/**
 * @author ecc4
 * <p>
 * To change the template for this generated type comment go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
public class IntegratorDCVGCMD extends IntegratorBox {

	IntegratorMC integratormc;
	IntegratorMD integratormd;
	public double zFraction = 0.1;
	private MyMCMove mcMove1, mcMove2, mcMove3, mcMove4;
	private ISpecies speciesA, speciesB;
	private final PotentialMasterHybrid potentialMasterHybrid;
	private int MDStepCount, MDStepRepetitions;

	public IntegratorDCVGCMD(PotentialMaster parent, double temperature,
							 ISpecies species1, ISpecies species2, Box box) {
		super(parent, temperature, box);
		this.speciesA = species1;
		this.speciesB = species2;
		potentialMasterHybrid = (parent instanceof PotentialMasterHybrid)
				? (PotentialMasterHybrid) parent : null;
		setMDStepRepetitions(50);
	}

	public void setMDStepRepetitions(int interval) {
		MDStepRepetitions = interval;
		if (MDStepCount > interval || MDStepCount == 0) MDStepCount = interval;
	}

	protected void setup() {
		super.setup();
		potentialMasterHybrid.setUseNbrLists(false);
		integratormc.reset();
		potentialMasterHybrid.setUseNbrLists(true);
		integratormd.reset();
	}

	public void setTemperature(double t) {
		super.setTemperature(t);
		if (integratormc != null) {
			integratormc.setTemperature(t);
		}
		if (integratormd != null) {
			integratormd.setTemperature(t);
		}
	}

	public void setIsothermal(boolean b) {
		super.setIsothermal(b);
		integratormd.setIsothermal(b);
	}

	protected void doStepInternal() {
		if (potentialMasterHybrid != null) {
			potentialMasterHybrid.setUseNbrLists(MDStepCount > 0);
		}
		if (MDStepCount == 0) {
			MDStepCount = MDStepRepetitions;
			mcMove1.setupActiveAtoms();
			mcMove2.setupActiveAtoms();
			mcMove3.setupActiveAtoms();
			mcMove4.setupActiveAtoms();
			for (int i = 0; i < 50; i++) {
				integratormc.doStep();
			}
			IAtomList allAtoms = box.getLeafList();
			for (int i = 0; i < allAtoms.size(); i++) {
				if (allAtoms.get(i).getPosition().getX(2) < -40) {
					throw new RuntimeException(i + " " + allAtoms.get(i) + " " + allAtoms.get(i).getPosition());
				}
			}
			potentialMasterHybrid.setUseNbrLists(true);
			potentialMasterHybrid.getNeighborManager(box).reset();
			integratormd.reset();
		} else {
			MDStepCount--;
			if (MDStepCount % 50 == 0) {
				integratormd.reset();
			}
			IAtomList allAtoms = box.getLeafList();
			double Lz = box.getBoundary().getBoxSize().getX(2);
			for (int i = 0; i < allAtoms.size(); i++) {
				if (Math.abs(allAtoms.get(i).getPosition().getX(2)) > Lz / 2
						|| Math.abs(((IAtomKinetic) allAtoms.get(i)).getVelocity().getX(2)) > 100) {
					IMolecule m = allAtoms.get(i).getParentGroup();
					IAtomList atoms = m.getChildList();
					System.out.println("step " + stepCount);
					for (IAtom a : atoms) {
						System.out.println(a + " " + a.getPosition().getX(2) + " " + ((IAtomKinetic) a).getVelocity().getX(2));
					}
					throw new RuntimeException("oops");
				}
			}
			if (false && box.getNMolecules(speciesB) > 0) {
				IAtom a = box.getMoleculeList(speciesB).get(0).getChildList().get(1);
				System.out.println("step " + stepCount + " " + a.getPosition().getX(2));
			}
			integratormd.doStep();
			for (int i = 0; i < allAtoms.size(); i++) {
				if (Math.abs(allAtoms.get(i).getPosition().getX(2)) > Lz / 2
						|| Math.abs(((IAtomKinetic) allAtoms.get(i)).getVelocity().getX(2)) > 100) {
					IMolecule m = allAtoms.get(i).getParentGroup();
					IAtomList atoms = m.getChildList();
					System.out.println("step " + stepCount);
					for (IAtom a : atoms) {
						System.out.println(a + " " + a.getPosition().getX(2) + " " + ((IAtomKinetic) a).getVelocity().getX(2));
					}
					throw new RuntimeException("oops");
				}
			}
		}
	}

	public double getCurrentTime() {
		return integratormd.getCurrentTime();
	}

	public void setIntegrators(IntegratorMC intmc, IntegratorMD intmd, IRandom random) {
		integratormc = intmc;
		integratormd = intmd;
		integratormc.setTemperature(temperature);
		integratormd.setTemperature(temperature);
		mcMove1 = new MyMCMove(this, random, space, -zFraction);
		mcMove1.setMu(Kelvin.UNIT.toSim(-10000));
		mcMove2 = new MyMCMove(this, random, space, +zFraction);
		mcMove2.setMu(Kelvin.UNIT.toSim(-4419));
		MCMoveManager moveManager = integratormc.getMoveManager();
		moveManager.addMCMove(mcMove1);
		moveManager.addMCMove(mcMove2);
		mcMove1.setSpecies(speciesA);
		mcMove2.setSpecies(speciesA);
		mcMove3 = new MyMCMove(this, random, space, -zFraction);
		mcMove3.setMu(Kelvin.UNIT.toSim(-10000));
		mcMove4 = new MyMCMove(this, random, space, +zFraction);
		mcMove4.setMu(Kelvin.UNIT.toSim(-4419));
		moveManager.addMCMove(mcMove3);
		moveManager.addMCMove(mcMove4);
		mcMove3.setSpecies(speciesB);
		mcMove4.setSpecies(speciesB);
	}

	public void reset() {
		super.reset();
		potentialMasterHybrid.setUseNbrLists(false);
		integratormc.reset();
		potentialMasterHybrid.setUseNbrLists(true);
		integratormd.reset();
	}

	public MyMCMove[] mcMoves() {
		return new MyMCMove[]{mcMove1, mcMove2, mcMove3, mcMove4};
	}
}
