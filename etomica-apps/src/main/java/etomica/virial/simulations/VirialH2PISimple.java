/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.virial.simulations;

import com.fasterxml.jackson.databind.ObjectMapper;
import etomica.atom.*;
import etomica.atom.iterator.ANIntragroupExchange;
import etomica.atom.iterator.ApiIntergroupCoupled;
import etomica.chem.elements.Hydrogen;
import etomica.config.ConformationLinear;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCovariance;
import etomica.data.AccumulatorRatioAverageCovarianceFull;
import etomica.data.IData;
import etomica.data.histogram.HistogramSimple;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.integrator.mcmove.MCMove;
import etomica.math.DoubleRange;
import etomica.molecule.IMoleculeList;
import etomica.potential.*;
import etomica.potential.P1HydrogenMielke.P1HydrogenMielkeAtomic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesSpheresHetero;
import etomica.units.BohrRadius;
import etomica.units.Kelvin;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.cluster.Standard;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;

public class VirialH2PISimple {
	public static final double massH = Hydrogen.INSTANCE.getMass();

	public static void main(String[] args) {
		VirialH2PISimpleParam params = new VirialH2PISimpleParam();
		boolean isCommandLine = args.length > 0;
		if (isCommandLine) {
			ParseArgs parseArgs = new ParseArgs(params);
			parseArgs.parseArgs(args, true);
		}
		else {
			// default options - choose these before committing to CVS
			params.nBeads = 8;
			params.temperatureK = 500;
			params.numSteps = (long)1E6;

			// runtime options - make changes in these and not the default options above
			//			params.nBeads = 8;
			//			params.temperatureK = 100;
			//			params.numSteps = (long)1E5;
			//			params.cont = new contribution[]{contribution.B, contribution.B};
			//			params.potentialLevel = levelOptions.HINDE_PATKOWSKI;
			//			params.blOption = blOptions.VARIABLE;
		}

		final int nPoints = params.nPoints;
		final double temperatureK = params.temperatureK;
		final double temperature = Kelvin.UNIT.toSim(temperatureK);

		long steps = params.numSteps;
		double sigmaHSRef = params.sigmaHSRef;
		double refFrac = params.refFrac;
		final int nBeads = params.nBeads;
		final int beadFac = params.beadFac;

		final double[] HSB = new double[8];
		HSB[2] = Standard.B2HS(sigmaHSRef);
		HSB[3] = Standard.B3HS(sigmaHSRef);
		HSB[4] = Standard.B4HS(sigmaHSRef);
		HSB[5] = Standard.B5HS(sigmaHSRef);
		HSB[6] = Standard.B6HS(sigmaHSRef);
		HSB[7] = Standard.B7HS(sigmaHSRef);

		final levelOptions potentialLevel = params.potentialLevel;
		final blOptions blOption = params.blOption;
		final contribution[] cont = params.cont;
		System.out.println("Potential level = "+potentialLevel);
		System.out.println("bl Option = "+blOption);
		System.out.print("Contribution = ");
		for (int i=0; i<nPoints; i++) System.out.print(cont[i]+" ");
		System.out.println();
		System.out.println("Temperature: "+temperatureK+"K"+ " nBeads: "+nBeads);
		System.out.println("Reference diagram: B"+nPoints+" for hard spheres with diameter " + sigmaHSRef + " Angstroms");
		System.out.println("B"+nPoints+"HS: "+HSB[nPoints]);
		if (steps%1000 != 0) {
			throw new RuntimeException("steps should be a multiple of 1000");
		}
		System.out.println(steps+" steps (1000 IntegratorOverlap steps of "+(steps/1000)+")");

		Space space = Space3D.getInstance();
		// make ref and tar clusters
		MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
		ClusterWheatleyHS refCluster = new ClusterWheatleyHS(nPoints, fRef);
		refCluster.setTemperature(temperature);

		final P2HydrogenPatkowskiIso p2patIso = new P2HydrogenPatkowskiIso(space);
		final IPotentialAtomic p2pat = new P2HydrogenPatkowskiAtomic(space);
		final IPotentialAtomic p2hp = new P2HydrogenHindePatkowskiAtomic(space);

		IPotentialAtomic p2 = null;

		if (potentialLevel == levelOptions.ISO) p2 = p2patIso;
		if (potentialLevel == levelOptions.PATKOWSKI) p2 = p2pat;
		if (potentialLevel == levelOptions.HINDE_PATKOWSKI) p2 = p2hp;

		PotentialGroupPI pTarGroup = new PotentialGroupPI(beadFac);
		pTarGroup.addPotential(p2, new ApiIntergroupCoupled());
		MayerGeneral fTar = new MayerGeneral(pTarGroup) {
			@Override
			public double f(IMoleculeList pair, double r2, double beta) {
				return super.f(pair, r2, beta/nBeads);
			}
		};

		ClusterWheatleySoft tarCluster = new ClusterWheatleySoft(nPoints, fTar, 1e-12);
		tarCluster.setTemperature(temperature);

		// make species
		AtomTypeOriented atype = new AtomTypeOriented(Hydrogen.INSTANCE, space);
		SpeciesSpheresHetero speciesH2 = null;
		if (blOption == blOptions.FIXED_GROUND) {
			speciesH2 = new SpeciesSpheresHetero(space, new AtomTypeOriented[]{atype}) {
				@Override
				protected IAtom makeLeafAtom(AtomType leafType) {
					double bl = BohrRadius.UNIT.toSim(1.448736);
					return new AtomHydrogen(space, (AtomTypeOriented) leafType, bl);
				}
			};
		}
		else {
			speciesH2 = new SpeciesSpheresHetero(space, new AtomTypeOriented[]{atype}) {
				@Override
				protected IAtom makeLeafAtom(AtomType leafType) {
					double bl = AtomHydrogen.getAvgBondLength(temperatureK);
					return new AtomHydrogen(space, (AtomTypeOriented) leafType, bl);
				}
			};
		}
		speciesH2.setChildCount(new int [] {nBeads});
		speciesH2.setConformation(new ConformationLinear(space, 0));

		// make simulation
		final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, speciesH2, temperature, refCluster, tarCluster);
		PotentialGroup pIntra1 = sim.integrators[1].getPotentialMaster().makePotentialGroup(1);
		int[] seeds = sim.getRandomSeeds();
		System.out.println("Random seeds: "+ Arrays.toString(seeds));
		P1HydrogenMielkeAtomic p1 = new P1HydrogenMielkeAtomic(space);
		pIntra1.addPotential(p1, new ANIntragroupExchange(1, nBeads));
		//        We use ANIntragroupExchange here by purpose even though we may not be doing exchange
		sim.integrators[1].getPotentialMaster().addPotential(pIntra1,new ISpecies[]{sim.getSpecies(0)});
		//        sim.init();
		sim.integratorOS.setNumSubSteps(1000);
		steps /= 1000;

		//        add additional moves here, simulation already has translation and rotation moves
		//        rotation is a bit pointless when we can regrow the ring completely

		if (nBeads != 1) {
			sim.integrators[0].getMoveManager().removeMCMove(sim.mcMoveRotate[0]);
			sim.integrators[1].getMoveManager().removeMCMove(sim.mcMoveRotate[1]);
		}
		//        remove translation too
		//		sim.integrators[0].getMoveManager().removeMCMove(sim.mcMoveTranslate[0]);
		//		sim.integrators[1].getMoveManager().removeMCMove(sim.mcMoveTranslate[1]);

		//        ring regrow translation
		MCMoveClusterRingRegrow refTr = new MCMoveClusterRingRegrow(sim.getRandom(), space);
		double lambda = Constants.PLANCK_H/Math.sqrt(2*Math.PI*massH*temperature);
		refTr.setEnergyFactor(2*nBeads*Math.PI/(lambda*lambda));

		MCMoveClusterRingRegrow tarTr = new MCMoveClusterRingRegrow(sim.getRandom(), space);
		tarTr.setEnergyFactor(2*nBeads*Math.PI/(lambda*lambda));

		if (nBeads != 1) {
			sim.integrators[0].getMoveManager().addMCMove(refTr);
			sim.integrators[1].getMoveManager().addMCMove(tarTr);
		}
		//        ring regrow orientation
		MCMoveClusterRingRegrowOrientation refOr = new MCMoveClusterRingRegrowOrientation(sim.getRandom(), space, nBeads);
		MCMoveClusterRingRegrowOrientation tarOr = new MCMoveClusterRingRegrowOrientation(sim.getRandom(), space, nBeads);
		refOr.setStiffness(temperature, massH/2.0);
		tarOr.setStiffness(temperature, massH/2.0);
		boolean[] xcFlag = new boolean[nPoints];
		for (int i=0; i<nPoints; i++) {
			xcFlag[i] = cont[i] == contribution.XC;
		}
		refOr.setDoExchange(xcFlag);
		tarOr.setDoExchange(xcFlag);
		boolean fixedOrientation = false;
		if (nBeads != 1) {
			sim.integrators[0].getMoveManager().addMCMove(refOr);
			sim.integrators[1].getMoveManager().addMCMove(tarOr);
		}
		else {
			fixedOrientation = true;
		}
		MCMoveChangeBondLength cbl0 = new MCMoveChangeBondLength(sim.integrators[0].getPotentialMaster(), sim.getRandom(),space,temperature);
		MCMoveChangeBondLength cbl1 = new MCMoveChangeBondLength(sim.integrators[0].getPotentialMaster(), sim.getRandom(),space,temperature);

		if (blOption == blOptions.VARIABLE) {
			sim.integrators[0].getMoveManager().addMCMove(cbl0);
			sim.integrators[1].getMoveManager().addMCMove(cbl1);
			cbl0.setFixedOrientation(fixedOrientation);
			cbl1.setFixedOrientation(fixedOrientation);
			cbl0.setDoExchange(xcFlag);
			cbl1.setDoExchange(xcFlag);
			cbl0.setStiffness(massH,p1);
			cbl1.setStiffness(massH,p1);
		}
		System.out.println();
		String refFileName = params.refPrefFileName;
		if (refFileName.isEmpty()) {
			String tempString = ""+temperatureK;
			if (temperatureK == (int)temperatureK) {
				// temperature is an integer, use "200" instead of "200.0"
				tempString = ""+(int)temperatureK;
			}
			refFileName = "refpref"+nPoints;
			refFileName += "_2b";
			refFileName += "_"+tempString+"_PI";
		}
		long t1 = System.currentTimeMillis();
		// if using 3-body potential for B3, we must select initial configuration
		// that is not overlapping for any two molecules


		for (int b=0; b<nPoints; b++) {
			IMoleculeList molecules = sim.box[b].getMoleculeList();
			for (int m = 0; m<molecules.size(); m++) {
				if (xcFlag[m]) {
					IAtomList atoms = molecules.get(m).getChildList();
					for (int i = 0; i<atoms.size(); i++) {
						AtomHydrogen o = (AtomHydrogen)atoms.get(i);
						double cT = Math.cos((Math.PI*i)/atoms.size());
						double sT = Math.sin((Math.PI*i)/atoms.size());
						Vector vec = space.makeVector();
						vec.setX(0, cT);
						vec.setX(1, sT);
						o.getOrientation().setDirection(vec);
					}
				}
			}
		}


		sim.initRefPref(refFileName, steps/40);
		if (refFrac >= 0) {
			sim.integratorOS.setRefStepFraction(refFrac);
			sim.integratorOS.setAdjustStepFraction(false);
		}
		sim.equilibrate(refFileName, steps/10);
		System.out.println("equilibration finished");

		final HistogramSimple h1 = new HistogramSimple(500, new DoubleRange(0,Math.PI));
		IntegratorListener histListenerTarget = new IntegratorListener() {
			@Override
			public void integratorInitialized(IntegratorEvent e) {}
			@Override
			public void integratorStepStarted(IntegratorEvent e) {}
			@Override
			public void integratorStepFinished(IntegratorEvent e) {
				IAtomList atoms = sim.box[1].getLeafList();
				Vector a0 = ((IAtomOriented)atoms.get(0)).getOrientation().getDirection();
				Vector a1 = ((IAtomOriented)atoms.get(1)).getOrientation().getDirection();
				double angle = Math.acos(a0.dot(a1));
				h1.addValue(angle);
			}
		};
		sim.integrators[1].getEventManager().addListener(histListenerTarget);

		sim.integratorOS.setNumSubSteps((int)steps);
		sim.setAccumulatorBlockSize(steps);
		sim.integratorOS.setAggressiveAdjustStepFraction(true);
		sim.ai.setMaxSteps(1000);

		System.out.println("MC Move step sizes (ref)    "+sim.mcMoveTranslate[0].getStepSize());
		System.out.println("MC Move step sizes (target) "+sim.mcMoveTranslate[1].getStepSize());

		final double refIntegralF = HSB[nPoints];
		if (! isCommandLine) {
			IntegratorListener progressReport = new IntegratorListener() {
				@Override
				public void integratorInitialized(IntegratorEvent e) {}
				@Override
				public void integratorStepStarted(IntegratorEvent e) {}
				@Override
				public void integratorStepFinished(IntegratorEvent e) {
					if ((sim.integratorOS.getStepCount()*10) % sim.ai.getMaxSteps() != 0) return;
					System.out.print(sim.integratorOS.getStepCount()+" steps: ");
					double[] ratioAndError = sim.dvo.getAverageAndError();
					double ratio = ratioAndError[0];
					double error = ratioAndError[1];
					//                    System.out.println("abs average: "+ratio*refIntegralF+" error: "+error*refIntegralF);
					System.out.println("abs average: "+ratio*refIntegralF+" error: "+error*refIntegralF);
					if (ratio == 0 || Double.isNaN(ratio)) {
						throw new RuntimeException("oops");
					}
				}
			};
			sim.integratorOS.getEventManager().addListener(progressReport);
		}
		// this is where the simulation takes place
		sim.getController().actionPerformed();
		//end of simulation
		long t2 = System.currentTimeMillis();

		System.out.println(" ");
		System.out.println("final reference step fraction "+sim.integratorOS.getIdealRefStepFraction());
		System.out.println("actual reference step fraction "+sim.integratorOS.getRefStepFraction());
		System.out.println(" ");

		System.out.println("Reference system: ");
		List<MCMove> refMoves = sim.integrators[0].getMoveManager().getMCMoves();
		for (MCMove m : refMoves) {
			double acc = m.getTracker().acceptanceRatio();
			System.out.println(m.toString()+" acceptance ratio: "+acc);
		}
		System.out.println("Target system: ");
		List<MCMove> tarMoves = sim.integrators[1].getMoveManager().getMCMoves();
		for (MCMove m : tarMoves) {
			double acc = m.getTracker().acceptanceRatio();

			//            if (acc == 1) {
			//                throw new RuntimeException("something seems fishy");
			//            }
			System.out.println(m.toString()+" acceptance ratio: "+acc);
		}
		// Printing results here
		double[] ratioAndError = sim.dvo.getAverageAndError();
		double ratio = ratioAndError[0];
		double error = ratioAndError[1];
		double bn = ratio*HSB[nPoints];
		double bnError = error*Math.abs(HSB[nPoints]);

		System.out.println("ratio average: "+ratio+" error: "+error);
		System.out.println("abs average: "+bn+" error: "+bnError);
		DataGroup allYourBase = (DataGroup)sim.accumulators[0].getData();
		IData ratioData = allYourBase.getData(AccumulatorRatioAverageCovarianceFull.RATIO.index);
		IData ratioErrorData = allYourBase.getData(AccumulatorRatioAverageCovarianceFull.RATIO_ERROR.index);
		IData averageData = allYourBase.getData(AccumulatorAverage.AVERAGE.index);
		IData stdevData = allYourBase.getData(AccumulatorAverage.STANDARD_DEVIATION.index);
		IData errorData = allYourBase.getData(AccumulatorAverage.ERROR.index);
		IData correlationData = allYourBase.getData(AccumulatorAverage.BLOCK_CORRELATION.index);
		IData covarianceData = allYourBase.getData(AccumulatorAverageCovariance.BLOCK_COVARIANCE.index);
		double correlationCoef = covarianceData.getValue(1)/Math.sqrt(covarianceData.getValue(0)*covarianceData.getValue(3));
		correlationCoef = (Double.isNaN(correlationCoef) || Double.isInfinite(correlationCoef)) ? 0 : correlationCoef;
		double refAvg = averageData.getValue(0);
		double refOvAvg = averageData.getValue(1);
		System.out.print(String.format("reference ratio average: %20.15e error:  %10.5e  cor: %6.4f\n", ratioData.getValue(1), ratioErrorData.getValue(1), correlationCoef));
		System.out.print(String.format("reference average: %20.15e stdev: %9.4e error: %9.4e cor: %6.4f\n",
				averageData.getValue(0), stdevData.getValue(0), errorData.getValue(0), correlationData.getValue(0)));
		System.out.print(String.format("reference overlap average: %20.15e stdev: %9.4e error: %9.3e cor: %6.4f\n",
				averageData.getValue(1), stdevData.getValue(1), errorData.getValue(1), correlationData.getValue(1)));

		allYourBase = (DataGroup)sim.accumulators[1].getData();
		ratioData = allYourBase.getData(AccumulatorRatioAverageCovarianceFull.RATIO.index);
		ratioErrorData = allYourBase.getData(AccumulatorRatioAverageCovarianceFull.RATIO_ERROR.index);
		averageData = allYourBase.getData(AccumulatorAverage.AVERAGE.index);
		stdevData = allYourBase.getData(AccumulatorAverage.STANDARD_DEVIATION.index);
		errorData = allYourBase.getData(AccumulatorAverage.ERROR.index);
		correlationData = allYourBase.getData(AccumulatorAverage.BLOCK_CORRELATION.index);
		covarianceData = allYourBase.getData(AccumulatorAverageCovariance.BLOCK_COVARIANCE.index);
		int n = sim.numExtraTargetClusters;
		correlationCoef = covarianceData.getValue(n+1)/Math.sqrt(covarianceData.getValue(0)*covarianceData.getValue((n+2)*(n+2)-1));
		correlationCoef = (Double.isNaN(correlationCoef) || Double.isInfinite(correlationCoef)) ? 0 : correlationCoef;
		double tarAvg = averageData.getValue(0);
		double tarOvAvg = averageData.getValue(1);
		double tarCorr = correlationData.getValue(0);
		System.out.print(String.format("target ratio average: %20.15e  error: %10.5e  cor: %6.4f\n", ratioData.getValue(n+1), ratioErrorData.getValue(n+1), correlationCoef));
		System.out.print(String.format("target average: %20.15e stdev: %9.4e error: %9.4e cor: %6.4f\n",
				averageData.getValue(0), stdevData.getValue(0), errorData.getValue(0), correlationData.getValue(0)));
		System.out.print(String.format("target overlap average: %20.15e stdev: %9.4e error: %9.4e cor: %6.4f\n",
				averageData.getValue(n+1), stdevData.getValue(n+1), errorData.getValue(n+1), correlationData.getValue(n+1)));
		if (isCommandLine) {

			LinkedHashMap<String, Comparable> resultsMap = new LinkedHashMap<String, Comparable>();
			resultsMap.put("temperature", temperatureK);
			resultsMap.put("P", nBeads);
			resultsMap.put("bn", bn);
			resultsMap.put("bnError", bnError);
			resultsMap.put("refAvg", refAvg);
			resultsMap.put("refOvAvg", refOvAvg);
			resultsMap.put("tarAvg", tarAvg);
			resultsMap.put("tarOvAvg", tarOvAvg);
			resultsMap.put("tarCorr", tarCorr);
			for (MCMove m : tarMoves) {
				double acc = m.getTracker().acceptanceRatio();
				resultsMap.put(m.toString(), acc);
			}

			if ((t2-t1)/1000.0 > 24*3600) {
				resultsMap.put("time",(t2-t1)/(24*3600*1000.0));
				resultsMap.put("unit","days");
				System.out.println("time: "+(t2-t1)/(24*3600*1000.0)+" days");
			}
			else if ((t2-t1)/1000.0 > 3600) {
				resultsMap.put("time",(t2-t1)/(3600*1000.0));
				resultsMap.put("unit","hrs");
				System.out.println("time: "+(t2-t1)/(3600*1000.0)+" hrs");
			}
			else if ((t2-t1)/1000.0 > 60) {
				resultsMap.put("time",(t2-t1)/(60*1000.0));
				resultsMap.put("unit","mins");
				System.out.println("time: "+(t2-t1)/(60*1000.0)+" mins");
			}
			else {
				resultsMap.put("time",(t2-t1)/(1000.0));
				resultsMap.put("unit","secs");
				System.out.println("time: "+(t2-t1)/1000.0+" secs");
			}

			try {
				FileWriter jsonFile = new FileWriter(params.jsonOutputFileName);
				ObjectMapper om = new ObjectMapper();
				jsonFile.write(om.writeValueAsString(resultsMap));
				jsonFile.write("\n");
				jsonFile.close();
			} catch (IOException e) {
				throw new RuntimeException(e.getMessage());
			}
		}
		//        sim.printResults(HSB[nPoints]);

		if ((t2-t1)/1000.0 > 24*3600) {
			System.out.println("time: "+(t2-t1)/(24*3600*1000.0)+" days");
		}
		else if ((t2-t1)/1000.0 > 3600) {
			System.out.println("time: "+(t2-t1)/(3600*1000.0)+" hrs");
		}
		else if ((t2-t1)/1000.0 > 60) {
			System.out.println("time: "+(t2-t1)/(60*1000.0)+" mins");
		}
		else {
			System.out.println("time: "+(t2-t1)/1000.0+" secs");
		}
	}
	enum levelOptions {
		ISO, PATKOWSKI, HINDE_PATKOWSKI
	}
	enum blOptions {
		FIXED_GROUND, FIXED_TEMPAVG, VARIABLE
	}
	enum contribution {
		// Boltzmann - B, Exchange - XC;
		B, XC
	}
	/**
	 * Inner class for parameters
	 */
	public static class VirialH2PISimpleParam extends ParameterBase {
		public int nPoints = 2;
		public int nBeads = 8;
		public double temperatureK = 500.0;   // Kelvin
		public long numSteps = 1000000;
		public double refFrac = -1;
		public double sigmaHSRef = 3.0; // -1 means use equation for sigmaHSRef
		public int beadFac = 2;
		public levelOptions potentialLevel = levelOptions.PATKOWSKI;
		public blOptions blOption = blOptions.FIXED_GROUND;
		public contribution[] cont = new contribution[nPoints];
		public String jsonOutputFileName = "";
		public String refPrefFileName = "";
	}
}
