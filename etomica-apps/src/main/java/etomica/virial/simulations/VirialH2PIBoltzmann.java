/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtomList;
import etomica.atom.iterator.ApiIntergroupCoupled;
import etomica.chem.elements.Hydrogen;
import etomica.config.ConformationLinear;
import etomica.data.AccumulatorRatioAverageCovariance;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.graphics.*;
import etomica.integrator.IntegratorListenerAction;
import etomica.molecule.IMoleculeList;
import etomica.potential.P1HydrogenMielke.P2HydrogenMielkeAtomic;
import etomica.potential.PotentialGroup;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;
import etomica.units.Kelvin;
import etomica.units.Pixel;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.*;
import etomica.util.Constants;
import etomica.util.Constants.CompassDirection;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;

import javax.swing.*;
import java.awt.*;


public class VirialH2PIBoltzmann {

    public static void main(String[] args) {
        VirialHePIParam params = new VirialHePIParam();
        boolean isCommandline = args.length > 0;
        if (isCommandline) {
            ParseArgs parseArgs = new ParseArgs(params);
            parseArgs.parseArgs(args, true);
        }
        final double r0 = params.r00;
        final double kIntra = params.kIntra0;
        final int nPoints = params.nPoints;
        double temperature = params.temperature;
        final long steps = params.numSteps;
        final double aRef = params.aRef;
        final double aTarget = params.aTarget;
        if (aTarget > aRef) {
            throw new RuntimeException("aTarget should be less than aRef");
        }
        final double bRef = 1-aRef;
        final double bTarget = 1-aTarget;
        final int nSpheres = (params.nSpheres > -1) ? 2*params.nSpheres : 2*((int)(1200/temperature) + 7);

        Space space = Space3D.getInstance();

        PotentialGroup pTargetGroup = new PotentialGroup(2);
        System.out.println("H2 Path Integral ("+nSpheres+"-mer chains) B"+nPoints+" at "+temperature+"K");
        System.out.println("perturbing from a="+aRef+" to "+aTarget);
        temperature = Kelvin.UNIT.toSim(temperature);
//        P2HePCKLJS p2 = new P2HePCKLJS(space);
        P2HydrogenMielkeAtomic p2 = new P2HydrogenMielkeAtomic(space);
        double pU0 = -nSpheres*55095.0720797867;
        MayerENonGeneral eRef = new MayerENonGeneral(pTargetGroup, pU0) {
            public double f(IMoleculeList pair, double r2, double beta) {
            	double r = Math.sqrt(r2);
            	double uH = Math.exp(-beta*0.5*kIntra*(r - r0)*(r - r0));
//            	System.out.println("ref: "+uH);
                return (aRef + bRef*super.f(pair, r2, beta/nSpheres))*(aRef*uH + bRef);
            }
        };
        MayerENonGeneral eTarget = new MayerENonGeneral(pTargetGroup, pU0) {
            public double f(IMoleculeList pair, double r2, double beta) {
            	double r = Math.sqrt(r2);
            	double uH = Math.exp(-beta*0.5*kIntra*(r - r0)*(r - r0));
//            	System.out.println("tar: "+uH);
                return (aTarget + bTarget*super.f(pair, r2, beta/nSpheres))*(aTarget*uH + bTarget);

            }
        };

        ClusterSum refCluster = new ClusterSum(new ClusterBonds[]{new ClusterBonds(2, new int[][][]{{{0,1}}})}, new double[]{1.0}, new MayerFunction[]{eRef});
        ClusterSum targetCluster = new ClusterSum(new ClusterBonds[]{new ClusterBonds(2, new int[][][]{{{0,1}}})}, new double[]{1.0}, new MayerFunction[]{eTarget});
        ClusterWeight samplingCluster = ClusterWeightAbs.makeWeightCluster(refCluster);


        // the cluster's temperature determines the factor multiplied in the exponential (f=e-1)
        // we want 1/(P*kT)
        targetCluster.setTemperature(temperature);
        refCluster.setTemperature(temperature);

        System.out.println(steps+" steps");
        double h2Mass = 2*Hydrogen.INSTANCE.getMass();
        double lambda = Constants.PLANCK_H/Math.sqrt(2*Math.PI*h2Mass*temperature);
        double energyFac = nSpheres*Math.PI/(lambda*lambda);
        SpeciesGeneral species = new SpeciesBuilder(space)
                .addCount(AtomType.simple("He", h2Mass), nSpheres)
                .withConformation(new ConformationLinear(space, 0))
                .build();
        // the temperature here goes to the integrator, which uses it for the purpose of intramolecular interactions
        // we handle that manually below, so just set T=1 here
        final SimulationVirial sim = new SimulationVirial(space, species, 1.0, samplingCluster, refCluster, new ClusterAbstract[]{targetCluster});

//        sim.integrator.getMoveManager().removeMCMove(sim.mcMoveTranslate);
        sim.integrator.getMoveManager().removeMCMove(sim.mcMoveRotate);

//        MCMoveBox ring;
//        double energyFac = nSpheres*Math.PI/(lambda*lambda);
//        if (aRef < 1) {
//            MCMoveClusterRingPartialRegrow ringMove = new MCMoveClusterRingPartialRegrow(sim.integrator.getPotentialMaster(), sim.getRandom(), space, new int[][]{{0,1}});
//            ringMove.setEnergyFactor(energyFac);
//            int numRegrowBeads = nSpheres/3;
//            System.out.println("regrow "+numRegrowBeads+" beads");
//            ringMove.setNumBeads(numRegrowBeads);
//            ring = ringMove;
//        }
//        else {
//            MCMoveClusterRingRegrow ringMove = new MCMoveClusterRingRegrow(sim.getRandom(), space, new int[][]{{0,1}});
//            ringMove.setEnergyFactor(energyFac);
//            int numCBMCTrials = nSpheres/40+10;
//            System.out.println("regrow full ring");
//            ringMove.setNumTrial(numCBMCTrials);
//            ring = ringMove;
//        }
//        sim.integrator.getMoveManager().addMCMove(ring);
        MCMoveClusterRingRegrow ring0 = new MCMoveClusterRingRegrow(sim.getRandom(), space);

        ring0.setEnergyFactor(energyFac);


        sim.integrator.getMoveManager().addMCMove(ring0);


//        MCMoveClusterRingScale ringScale = new MCMoveClusterRingScale(sim.integrator.getPotentialMaster(), sim.getRandom(), space, new int[][]{{0,1}});
//        ringScale.setEnergyFactor(energyFac);
//        sim.integrator.getMoveManager().addMCMove(ringScale);
//        sim.integrator.getMoveManager().setFrequency(ringScale, 0.01);
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


        AtomType type = species.getAtomType(0);
        pTargetGroup.addPotential(p2, new ApiIntergroupCoupled());



        IAtomList leafList = sim.box.getLeafList();
        int half = leafList.size()/2;
        for (int i = half; i<leafList.size(); i++) {
            Vector p = leafList.get(i).getPosition();
            p.setX(2, r0);
        }
        sim.box.trialNotify();
        sim.box.acceptNotify();

        if (false) {
            double vSize = 10;
            sim.box.getBoundary().setBoxSize(Vector.of(new double[]{vSize, vSize, vSize}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
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
            ColorSchemeRandomByMolecule colorScheme = new ColorSchemeRandomByMolecule(sim.getSpeciesManager(), sim.box, sim.getRandom());
            displayBox.setColorScheme(colorScheme);
            simGraphic.makeAndDisplayFrame();


            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file

            final DisplayTextBox averageBox = new DisplayTextBox();
            averageBox.setLabel("Average");
            final DisplayTextBox errorBox = new DisplayTextBox();
            errorBox.setLabel("Error");
            JLabel jLabelPanelParentGroup = new JLabel("ratio");
            final JPanel panelParentGroup = new JPanel(new BorderLayout());
            panelParentGroup.add(jLabelPanelParentGroup,CompassDirection.NORTH.toString());
            panelParentGroup.add(averageBox.graphic(), BorderLayout.WEST);
            panelParentGroup.add(errorBox.graphic(), BorderLayout.EAST);
            simGraphic.getPanel().controlPanel.add(panelParentGroup, SimulationPanel.getVertGBC());


            IAction pushAnswer = new IAction() {
                DataDouble data = new DataDouble();

                public void actionPerformed() {
                    DataGroup allYourBase = (DataGroup)sim.accumulator.getData();
                    data.x = ((DataDoubleArray) allYourBase.getData(AccumulatorRatioAverageCovariance.RATIO.index)).getData()[1];
                    averageBox.putData(data);
                    data.x = ((DataDoubleArray) allYourBase.getData(AccumulatorRatioAverageCovariance.RATIO_ERROR.index)).getData()[1];
                    errorBox.putData(data);
                }
            };
            IDataInfo dataInfo = new DataDouble.DataInfoDouble("B"+nPoints, new CompoundDimension(new Dimension[]{new DimensionRatio(Volume.DIMENSION, Quantity.DIMENSION)}, new double[]{nPoints-1}));
            averageBox.putDataInfo(dataInfo);
            averageBox.setLabel("average");
            errorBox.putDataInfo(dataInfo);
            errorBox.setLabel("error");
            errorBox.setPrecision(2);
            sim.integrator.getEventManager().addListener(new IntegratorListenerAction(pushAnswer));

            sim.addEquilibration(steps / 100);
            sim.getController().addActivity(new ActivityIntegrate(sim.integrator));

            return;
        }
        
        
        sim.equilibrate(steps/100);
ActivityIntegrate ai = new ActivityIntegrate(sim.integrator, steps);
sim.setAccumulatorBlockSize(steps > 1000 ? steps/1000 : 1);

        System.out.println("equilibration finished");

//        if (false) {
//            final double refIntegralF = refIntegral;
//            IntegratorListener progressReport = new IntegratorListener() {
//                public void integratorInitialized(IIntegratorEvent e) {}
//                public void integratorStepStarted(IIntegratorEvent e) {}
//                public void integratorStepFinished(IIntegratorEvent e) {
//                    if ((sim.integrator.getStepCount()*10) % sim.getController2().getMaxSteps() != 0) return;
//                    System.out.print(sim.integrator.getStepCount()+" steps: ");
//                    double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
//                    double ratio = ratioAndError[0];
//                    double error = ratioAndError[1];
//                    System.out.println("abs average: "+ratio*refIntegralF+", error: "+error*refIntegralF);
//                }
//            };
//            sim.integratorOS.getEventManager().addListener(progressReport);
//        }

        sim.integrator.getMoveManager().setEquilibrating(false);
sim.getController().runActivityBlocking(ai);


        System.out.println("Ring acceptance "+ring0.getTracker().acceptanceRatio());
        if (aRef == 1) {
            double refIntegral = Math.pow(lambda, 3.0) * Math.pow(2.0, -2.5);
            
            System.out.println("reference integral "+refIntegral);
        }

        DataGroup allYourBase = (DataGroup)sim.accumulator.getData();
        
        System.out.println();
        System.out.println("reference average: " + ((DataDoubleArray) allYourBase.getData(sim.accumulator.AVERAGE.index)).getData()[0]
                + " stdev: " + ((DataDoubleArray) allYourBase.getData(sim.accumulator.STANDARD_DEVIATION.index)).getData()[0]
                + " error: " + ((DataDoubleArray) allYourBase.getData(sim.accumulator.ERROR.index)).getData()[0]);

        double ratio = ((DataDoubleArray) allYourBase.getData(sim.accumulator.RATIO.index)).getData()[1];
        double error = ((DataDoubleArray) allYourBase.getData(sim.accumulator.RATIO_ERROR.index)).getData()[1];

        System.out.println("target average: " + ((DataDoubleArray) allYourBase.getData(sim.accumulator.AVERAGE.index)).getData()[1]
                + " stdev: " + ((DataDoubleArray) allYourBase.getData(sim.accumulator.STANDARD_DEVIATION.index)).getData()[1]
                + " error: " + ((DataDoubleArray) allYourBase.getData(sim.accumulator.ERROR.index)).getData()[1]);

        System.out.println();
        System.out.println("ratio average: "+ratio+", error: "+error);
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

    /**
     * Inner class for parameters
     */
    public static class VirialHePIParam extends ParameterBase {
        public int nPoints = 2;
        public int nSpheres = -1;
        public double temperature = 500;   // Kelvin
        public long numSteps = 1000;
        public double aRef = 0.001;       // a=1 => ideal gas,  a=0 => ebond
        public double aTarget = 0;
        protected double r00 = 0.7414197982501762, kIntra0 = 346924.86139634805;
    }
}
