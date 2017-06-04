/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import java.awt.Color;
import java.util.List;

import javax.swing.JLabel;
import javax.swing.JPanel;

import etomica.action.IAction;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.IAtomType;
import etomica.api.IIntegratorEvent;
import etomica.api.IIntegratorListener;
import etomica.space.Vector;
import etomica.atom.AtomHydrogen;
import etomica.atom.AtomTypeOrientedSphere;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtomTypeOriented;
import etomica.chem.elements.Hydrogen;
import etomica.config.ConformationLinear;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.IEtomicaDataInfo;
import etomica.data.meter.MeterAvgBondLength;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorSchemeRandomByMolecule;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.integrator.mcmove.MCMove;
import etomica.listener.IntegratorListenerAction;
import etomica.potential.P1HydrogenMielke.P1HydrogenMielkeAtomic;
import etomica.potential.P1IntraMolecular;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresHetero;
import etomica.units.BohrRadius;
import etomica.units.CompoundDimension;
import etomica.units.Dimension;
import etomica.units.DimensionRatio;
import etomica.units.Kelvin;
import etomica.units.Pixel;
import etomica.units.Quantity;
import etomica.units.Volume;
import etomica.util.Constants;
import etomica.util.Constants.CompassDirection;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterBonds;
import etomica.virial.ClusterExchange;
import etomica.virial.ClusterWeight;
import etomica.virial.ClusterWeightAbs;
import etomica.virial.MCMoveChangeBondLength;
import etomica.virial.MCMoveClusterRingRegrow;
import etomica.virial.MCMoveClusterRingRegrowExchange;
import etomica.virial.MCMoveClusterRingRegrowOrientation;
import etomica.virial.P1IntraLambda;


public class VirialH2PIXC {

	public static void main(String[] args) {
		VirialHePIParam params = new VirialHePIParam();
		boolean isCommandline = args.length > 0;
		if (isCommandline) {
			ParseArgs parseArgs = new ParseArgs(params);
			parseArgs.parseArgs(args, true);
		}
		final int nPoints = params.nPoints;
		final double temperatureK = params.temperature;
		final long steps = params.numSteps;
		final double aRef = params.aRef;
		final double aTarget = params.aTarget;
		if (aTarget < aRef) {
			throw new RuntimeException("aTarget should be more than aRef");
		}
		final double bRef = 1-aRef;
		final double bTarget = 1-aTarget;
		int nSpheres = (params.nSpheres > -1) ? 2*params.nSpheres : 2*((int)(1200/temperatureK) + 7);
		int exponent = (int)(Math.log(nSpheres)/Math.log(2)) + 1;
		nSpheres = (int)Math.pow(2, exponent);
		//        nSpheres = 64;

		Space space = Space3D.getInstance();


		System.out.println("H2 Path Integral ("+nSpheres+"-mer chains) B"+nPoints+" at "+temperatureK+"K");
		System.out.println("perturbing from a="+aRef+" to "+aTarget);
		final double temperature = Kelvin.UNIT.toSim(temperatureK);
		//        P2HePCKLJS p2 = new P2HePCKLJS(space);

		P1HydrogenMielkeAtomic p1 = new P1HydrogenMielkeAtomic(space);
		double r0 = 0.7414757777367901;
		double u0 = p1.u(r0);
		P1IntraLambda modifiedRefPot = new P1IntraLambda(space, aRef, p1,u0);
		P1IntraLambda modifiedTarPot = new P1IntraLambda(space, aTarget, p1,u0);

		//        MayerENonGeneral eRef = new MayerENonGeneral(pTargetGroup, pU0) {
		//            public double f(IMoleculeList pair, double r2, double beta) {
		//                return aRef + bRef*super.f(pair, r2, beta/nSpheres);
		//            }
		//        };
		//        MayerENonGeneral eTarget = new MayerENonGeneral(pTargetGroup, pU0) {
		//            public double f(IMoleculeList pair, double r2, double beta) {
		//                return aTarget + bTarget*super.f(pair, r2, beta/nSpheres);
		//            }
		//        };

		//        ClusterSum refCluster = new ClusterSum(new ClusterBonds[]{new ClusterBonds(2, new int[][][]{{{0,1}}})}, new double[]{1.0}, new MayerFunction[]{eRef});
		//        ClusterSum targetCluster = new ClusterSum(new ClusterBonds[]{new ClusterBonds(2, new int[][][]{{{0,1}}})}, new double[]{1.0}, new MayerFunction[]{eTarget});
		ClusterExchange refCluster = new ClusterExchange(modifiedRefPot);
		ClusterExchange targetCluster = new ClusterExchange(modifiedTarPot);
		ClusterWeight samplingCluster = ClusterWeightAbs.makeWeightCluster(refCluster);


		// the cluster's temperature determines the factor multiplied in the exponential (f=e-1)
		// we want 1/(P*kT)
		targetCluster.setTemperature(temperature);
		refCluster.setTemperature(temperature);

		System.out.println(steps+" steps");
		double hMass = Hydrogen.INSTANCE.getMass();
		//        System.out.println(hMass);
		//        System.exit(1);
		double lambda = Constants.PLANCK_H/Math.sqrt(2*Math.PI*hMass*temperature);
		//        SpeciesSpheres species = new SpeciesSpheres(space, nSpheres, new AtomTypeLeaf(new ElementChemical("He", h2Mass, 2)), new ConformationLinear(space, 0));
		IAtomTypeOriented atype = new AtomTypeOrientedSphere(Hydrogen.INSTANCE, space);
		SpeciesSpheresHetero species = new SpeciesSpheresHetero(space,new IAtomTypeOriented [] {atype}) {
			@Override
			protected IAtom makeLeafAtom(IAtomType leafType) {
				double bl = BohrRadius.UNIT.toSim(1.448736);// AtomHydrogen.getAvgBondLength(temperatureK);
				return new AtomHydrogen(space,(IAtomTypeOriented)leafType,bl);
			}
		};
		species.setChildCount(new int [] {nSpheres});
		species.setConformation(new ConformationLinear(space, 0));

		// the temperature here goes to the integrator, which uses it for the purpose of intramolecular interactions
		// we handle that manually below, so just set T=1 here
		//        int[] seed = new int [1];
		//        seed[0] = 1;
		final SimulationVirial sim = new SimulationVirial(space, species, 1.0, samplingCluster, refCluster, new ClusterAbstract[]{targetCluster});

		sim.integrator.getMoveManager().removeMCMove(sim.mcMoveTranslate);
		sim.integrator.getMoveManager().removeMCMove(sim.mcMoveRotate);

		double kHarmonic = nSpheres*Math.PI/(lambda*lambda);
		double energyFac = 2*kHarmonic;
		boolean[] flag = new boolean[nPoints];
		for (int i=0; i<nPoints; i++) flag[i] = true;
		MCMoveClusterRingRegrow ring0 = new MCMoveClusterRingRegrow(sim.getRandom(), space);
		ring0.setEnergyFactor(energyFac);
		MCMoveClusterRingRegrowOrientation move0 = new MCMoveClusterRingRegrowOrientation(sim.getRandom(), space, nSpheres);
		move0.setStiffness(temperature, species.getAtomType(0).getMass());

		move0.setDoExchange(flag);
		final MCMoveChangeBondLength cbl0 = new MCMoveChangeBondLength(sim.integrator.getPotentialMaster(), sim.getRandom(),space,temperature);
		cbl0.setDoExchange(flag);
		MCMoveClusterRingRegrowExchange exMove0 = new MCMoveClusterRingRegrowExchange(sim.getRandom(), space);
		exMove0.setEnergyFactor(kHarmonic);

		//        MCMoveChangeBondLengthBruteForce cblBF0 = new MCMoveChangeBondLengthBruteForce(sim.integrator.getPotentialMaster(), sim.getRandom(),space,temperature);
		//        cblBF0.setFixedOrientation(true);
		//        cblBF0.setDoExchange(false);

		sim.integrator.getMoveManager().addMCMove(ring0);
		sim.integrator.getMoveManager().addMCMove(move0);
		sim.integrator.getMoveManager().addMCMove(cbl0);
		cbl0.setStiffness(Hydrogen.INSTANCE.getMass(), modifiedRefPot);
		//        sim.integrator.getMoveManager().addMCMove(exMove0);
		//        sim.integrator.getMoveManager().addMCMove(cblBF0);
		//        ((MCMoveStepTracker)cblBF0.getTracker()).setNoisyAdjustment(true);
		//        cblBF0.setStepSize(0.05);
		//        cblBF0.setStiffness(species.getAtomType(0));


		//        AtomActionTranslateBy translator = new AtomActionTranslateBy(space);
		//        IVectorRandom groupTranslationVector = (IVectorRandom)translator.getTranslationVector();
		//        MoleculeChildAtomAction moveMoleculeAction = new MoleculeChildAtomAction(translator);
		//        IMoleculeList molecules = sim.box.getMoleculeList();
		//        for (int i=1; i<nPoints; i++) {
		//            groupTranslationVector.setX(0, i*3);
		//            moveMoleculeAction.actionPerformed(molecules.getMolecule(i));
		//            molecules.getMolecule(i).getChildList().getAtom(i).getPosition().setX(1, 1);
		//        }
		//        sim.box.trialNotify();
		//        double pi = sim.box.getSampleCluster().value(sim.box);
		//        if (pi == 0) throw new RuntimeException("initialization failed");
		//        sim.box.acceptNotify();

		IAtomType type = species.getAtomType(0);

		IAtomList leafList = sim.box.getLeafList();
		for (int i=0; i<leafList.getAtomCount(); i++) {

			AtomHydrogen o = (AtomHydrogen)leafList.getAtom(i);
			double cT = Math.cos((Math.PI*i)/leafList.getAtomCount());
			double sT = Math.sin((Math.PI*i)/leafList.getAtomCount());
			Vector vec = space.makeVector();
			vec.setX(0, cT);
			vec.setX(1, sT);
			o.getOrientation().setDirection(vec);
		}
		double newr0 = newtonRaphson(temperature, nSpheres, kHarmonic, Math.cos(Math.PI/nSpheres),p1);
		//        System.out.println(newr0);
		//        r0 = 0.7414757777367901;
		//        System.exit(1);

		for (int i=0; i<nSpheres; i++) {
			((AtomHydrogen)leafList.getAtom(i)).setBondLength(newr0);
		}
		sim.box.trialNotify();
		sim.box.acceptNotify();
		//        System.out.println(" "+sim.box.getSampleCluster().value(sim.box));

		if (false) {
			double vSize = 10;
			sim.box.getBoundary().setBoxSize(space.makeVector(new double[]{vSize,vSize,vSize}));
			SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space, sim.getController());
			DisplayBox displayBox = simGraphic.getDisplayBox(sim.box);
			displayBox.setPixelUnit(new Pixel(300.0/vSize));
			displayBox.setShowBoundary(false);
			((DisplayBoxCanvasG3DSys)displayBox.canvas).setBackgroundColor(Color.WHITE);

			//            IAtomList leafList = sim.box.getLeafList();
			//            AtomPair pair = new AtomPair();
			//            for (int i=0; i<leafList.getAtomCount()-1; i++) {
			//                pair.atom0 = leafList.getAtom(i);
			//                pair.atom1 = leafList.getAtom(i+1);
			//                ((DisplayBoxCanvasG3DSys)displayBox.canvas).makeBond(pair, null);
			//            }
			//            pair.atom0 = leafList.getAtom(leafList.getAtomCount()-1);
			//            pair.atom1 = leafList.getAtom(0);
			//            ((DisplayBoxCanvasG3DSys)displayBox.canvas).makeBond(pair, null);

			DiameterHashByType diameterManager = (DiameterHashByType)displayBox.getDiameterHash();
			diameterManager.setDiameter(type, 1.0/nSpheres);
			ColorSchemeRandomByMolecule colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box, sim.getRandom());
			displayBox.setColorScheme(colorScheme);
			simGraphic.makeAndDisplayFrame();


			// if running interactively, set filename to null so that it doens't read
			// (or write) to a refpref file

			final DisplayTextBox averageBox = new DisplayTextBox();
			averageBox.setLabel("Average");
			final DisplayTextBox errorBox = new DisplayTextBox();
			errorBox.setLabel("Error");
			JLabel jLabelPanelParentGroup = new JLabel("ratio");
			final JPanel panelParentGroup = new JPanel(new java.awt.BorderLayout());
			panelParentGroup.add(jLabelPanelParentGroup,CompassDirection.NORTH.toString());
			panelParentGroup.add(averageBox.graphic(), java.awt.BorderLayout.WEST);
			panelParentGroup.add(errorBox.graphic(), java.awt.BorderLayout.EAST);
			simGraphic.getPanel().controlPanel.add(panelParentGroup, SimulationPanel.getVertGBC());


			IAction pushAnswer = new IAction() {
				@Override
				public void actionPerformed() {
					DataGroup allYourBase = (DataGroup)sim.accumulator.getData();
					data.x = ((DataDoubleArray)allYourBase.getData(sim.accumulator.RATIO.index)).getData()[1];;
					averageBox.putData(data);
					data.x = ((DataDoubleArray)allYourBase.getData(sim.accumulator.RATIO_ERROR.index)).getData()[1];
					errorBox.putData(data);
				}

				DataDouble data = new DataDouble();
			};
			IEtomicaDataInfo dataInfo = new DataDouble.DataInfoDouble("B"+nPoints, new CompoundDimension(new Dimension[]{new DimensionRatio(Volume.DIMENSION, Quantity.DIMENSION)}, new double[]{nPoints-1}));
			averageBox.putDataInfo(dataInfo);
			averageBox.setLabel("average");
			errorBox.putDataInfo(dataInfo);
			errorBox.setLabel("error");
			errorBox.setPrecision(2);
			sim.integrator.getEventManager().addListener(new IntegratorListenerAction(pushAnswer));

			sim.getController().removeAction(sim.ai);
			sim.getController().addAction(new IAction() {
				@Override
				public void actionPerformed() {
					sim.equilibrate(steps/100);
					sim.ai.setMaxSteps(Long.MAX_VALUE);
				}
			});
			sim.getController().addAction(sim.ai);

			return;
		}


		sim.equilibrate(steps/100);

		sim.setAccumulatorBlockSize(steps > 1000 ? steps/1000 : 1);

		System.out.println("equilibration finished");
		final AccumulatorAverageFixed avg0 = new AccumulatorAverageFixed(steps);
		final MeterAvgBondLength mbl = new MeterAvgBondLength();
		mbl.setBox(sim.box);
		DataPumpListener pumpDiagrams1 = new DataPumpListener(mbl, avg0);
		sim.integrator.getEventManager().addListener(pumpDiagrams1);

		if (true) {
			IIntegratorListener progressReport = new IIntegratorListener() {
				@Override
				public void integratorInitialized(IIntegratorEvent e) {}
				@Override
				public void integratorStepStarted(IIntegratorEvent e) {}
				@Override
				public void integratorStepFinished(IIntegratorEvent e) {
					if ((sim.integrator.getStepCount()*10) % sim.ai.getMaxSteps() != 0) return;
					System.out.println(temperatureK+" "+avg0.getData(avg0.AVERAGE).getValue(0)+" "+avg0.getData(avg0.ERROR).getValue(0)+" "+avg0.getData(avg0.BLOCK_CORRELATION).getValue(0));
				}
			};
			sim.integrator.getEventManager().addListener(progressReport);
		}

		sim.integrator.getMoveManager().setEquilibrating(false);
		sim.ai.setMaxSteps(steps);
		sim.getController().actionPerformed();

		if (aRef == 1) {
			double refIntegral = Math.pow(lambda, 3.0) * Math.pow(2.0, -2.5);

			System.out.println("reference integral "+refIntegral);
		}

		DataGroup allYourBase = (DataGroup)sim.accumulator.getData();
		List<MCMove> moves = sim.integrator.getMoveManager().getMCMoves();
		for (MCMove m : moves) {
			System.out.println(m.toString()+" acceptance ratio : "+m.getTracker().acceptanceRatio());
		}
		//        System.out.println("Target Ring acceptance "+ring0.getTracker().acceptanceRatio());
		//        System.out.println("Orientation acceptance "+move0.getTracker().acceptanceRatio());
		//        System.out.println("Bond length acceptance "+cbl0.getTracker().acceptanceRatio());

		System.out.println();
		System.out.println("reference average: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.AVERAGE.index)).getData()[0]
				+" stdev: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.STANDARD_DEVIATION.index)).getData()[0]
						+" error: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.ERROR.index)).getData()[0]);

		double ratio = ((DataDoubleArray)allYourBase.getData(sim.accumulator.RATIO.index)).getData()[1];
		double error = ((DataDoubleArray)allYourBase.getData(sim.accumulator.RATIO_ERROR.index)).getData()[1];

		System.out.println("target average: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.AVERAGE.index)).getData()[1]
				+" stdev: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.STANDARD_DEVIATION.index)).getData()[1]
						+" error: "+((DataDoubleArray)allYourBase.getData(sim.accumulator.ERROR.index)).getData()[1]);

		System.out.println();
		System.out.println("ratio average: "+ratio+", error: "+error);

		//        for (int i=0; i<nSpheres; i++) {
		//        	double xAvg = cbl0.eAvg[i]/(cbl0.moveCount - 1000.0);
		//        	double x2Avg = cbl0.e2Avg[i]/(cbl0.moveCount - 1000.0);
		//        	double sd2 = x2Avg - xAvg*xAvg;
		//        	System.out.println(i+" actual = "+sd2+" generated = "+cbl0.sigma[i]*cbl0.sigma[i]);
		//        }
	}

	public static ClusterBonds[] append(ClusterBonds[] inArray, ClusterBonds[] newBonds) {
		ClusterBonds[] outArray = new ClusterBonds[inArray.length + newBonds.length];
		System.arraycopy(inArray, 0, outArray, 0, inArray.length);
		System.arraycopy(newBonds, 0, outArray, inArray.length, newBonds.length);
		return outArray;
	}

	public static double[] append(double[] inArray, double[] newWeights) {
		double[] outArray = new double[inArray.length + newWeights.length];
		System.arraycopy(inArray, 0, outArray, 0, inArray.length);
		System.arraycopy(newWeights, 0, outArray, inArray.length, newWeights.length);
		return outArray;
	}

	protected static double newtonRaphson(double t, int P, double kHarmonic, double cT, P1IntraMolecular p1) {
		int maxIter = 100000;
		double tol = 1E-15;
		double x0 = BohrRadius.UNIT.toSim(1.0);
		double xNew = 0;
		boolean done = false;
		double r0 = 0;
		for (int i=0; i<maxIter && !done; i++) {
			double d2vdr2 = 2*kHarmonic + 2/(P*x0*x0) + p1.d2u(x0)/(x0*x0*t*P);
			if (d2vdr2 == 0 || d2vdr2 != d2vdr2) throw new RuntimeException("x0 = "+x0+" d2vdr2 = "+d2vdr2);
			double dvdr = 2*x0*kHarmonic*(1 - cT) - 2.0/(P*x0) + p1.du(x0)/(t*P) ;
			xNew = x0 - dvdr/d2vdr2;
			if (xNew != xNew) throw new RuntimeException("xNew =" + xNew+ " x0 = "+x0+" dvdr = "+dvdr+" d2vdr2 ="+d2vdr2);
			if (xNew < 0) xNew = x0;
			if ((x0 - xNew)*(x0 - xNew) < tol*tol) {
				r0 = xNew;
				done = true;
			}
			x0 = xNew;
		}
		if (r0 == 0 || Double.isInfinite(r0) || Double.isNaN(r0)) throw new RuntimeException("Newton Raphson failed! "+ r0);
		return r0;
	}

	/**
	 * Inner class for parameters
	 */
	public static class VirialHePIParam extends ParameterBase {
		public int nPoints = 2;
		public int nSpheres = -1;
		public double temperature = 15;   // Kelvin
		public long numSteps = 1000000;
		public double aRef = 0.8;
		public double aTarget = 0.9;
	}
}

