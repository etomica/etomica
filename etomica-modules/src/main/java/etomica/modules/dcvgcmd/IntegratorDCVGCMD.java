/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/*
 * Created on Apr 12, 2004
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package etomica.modules.dcvgcmd;

import etomica.box.Box;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorMD;
import etomica.integrator.mcmove.MCMoveManager;
import etomica.modifier.Modifier;
import etomica.potential.compute.PotentialCompute;
import etomica.species.ISpecies;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Null;
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
	double zFraction = 0.1;
	private MyMCMove mcMove1, mcMove2, mcMove3, mcMove4;
	private final ISpecies speciesA, speciesB;
	private int MDStepCount, MDStepRepetitions;

	public IntegratorDCVGCMD(PotentialCompute compute, double temperature,
							 ISpecies species1, ISpecies species2, Box box) {
		super(compute, temperature, box);
		this.speciesA = species1;
		this.speciesB = species2;
		setMDStepRepetitions(50);
	}

	public void setMDStepRepetitions(int interval) {
		MDStepRepetitions = interval;
		if (MDStepCount > interval || MDStepCount == 0) MDStepCount = interval;
	}

	protected void setup() {
		super.setup();
		integratormc.reset();
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
		if (MDStepCount == 0) {
			MDStepCount = MDStepRepetitions;
			mcMove1.setupActiveAtoms();
			mcMove2.setupActiveAtoms();
			mcMove3.setupActiveAtoms();
			mcMove4.setupActiveAtoms();
			integratormc.reset();
			for (int i = 0; i < 50; i++) {
				integratormc.doStep();
			}
			integratormd.reset();
		} else {
			MDStepCount--;
			integratormd.doStep();
		}
	}

	//modulator for Mu's
	public class Mu1Modulator implements Modifier {
		public void setValue(double x) {
			setMu(x, mcMove2.getMu());
		}

		public String getLabel() {
			return "mu1";
		}

		public double getValue() {
			return mcMove1.getMu();
		}

		public Dimension getDimension() {
			return Null.DIMENSION;
		}
	}

	public class Mu2Modulator implements Modifier {
		public void setValue(double x) {
			setMu(mcMove1.getMu(), x);
		}

		public double getValue() {
			return mcMove2.getMu();
		}

		public Dimension getDimension() {
			return Null.DIMENSION;
		}

		public String getLabel() {
			return "mu2";
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
		mcMove2 = new MyMCMove(this, random, space, +zFraction);
		MCMoveManager moveManager = integratormc.getMoveManager();
		moveManager.addMCMove(mcMove1);
		moveManager.addMCMove(mcMove2);
		mcMove1.setSpecies(speciesA);
		mcMove2.setSpecies(speciesA);
		mcMove3 = new MyMCMove(this, random, space, -zFraction);
		mcMove4 = new MyMCMove(this, random, space, +zFraction);
		moveManager.addMCMove(mcMove3);
		moveManager.addMCMove(mcMove4);
		mcMove3.setSpecies(speciesB);
		mcMove4.setSpecies(speciesB);
	}

	public void setMu(double mu1, double mu2) {
		mcMove1.setMu(mu1);
		mcMove2.setMu(mu2);
		mcMove3.setMu(mu2);
		mcMove4.setMu(mu1);
	}

	public void reset() {
		super.reset();
		integratormc.reset();
		integratormd.reset();
	}

	public MyMCMove[] mcMoves() {
		return new MyMCMove[]{mcMove1, mcMove2, mcMove3, mcMove4};
	}
}
