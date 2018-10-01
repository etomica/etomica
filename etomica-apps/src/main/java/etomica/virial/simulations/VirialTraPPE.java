/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.IAction;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.iterator.ApiBuilder;
import etomica.box.Box;
import etomica.chem.elements.*;
import etomica.config.ConformationGeneric;
import etomica.config.IConformation;
import etomica.data.histogram.HistogramSimple;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.math.DoubleRange;
import etomica.math.SpecialFunctions;
import etomica.molecule.IMoleculeList;
import etomica.potential.*;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresHetero;
import etomica.units.Degree;
import etomica.units.Electron;
import etomica.units.Kelvin;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;
import etomica.virial.*;
import etomica.virial.cluster.Standard;

import java.awt.*;
import java.util.Arrays;

/**
 * Mayer sampling simulation for alkanes using the TraPPE force field.
 *   M.G. Martin and J.I. Siepmann, "Transferable Potentials for Phase
 *   Equilibria. 1. United-Atom Description of n-Alkanes," J. Phys. Chem. B
 *   102, 2569-2577 (1998)
 */
public class VirialTraPPE {


    public static void main(String[] args) {
        VirialParam params = new VirialParam();
        boolean isCommandline = args.length > 0;
        if (isCommandline) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            // customize parameters here
            params.chemForm = ChemForm.N2;
            params.nPoints = 2;
            params.nDer = 0;
            params.temperature = 400;
            params.numSteps = 1000000;

            params.refFrac = -1;
            params.sigmaHSRef = 5;
            params.seed = null;

            params.doHist = false;
            params.dorefpref = false;
            params.doChainRef = true;

            params.BDtol = 1e-12;
        }

        final ChemForm chemForm = params.chemForm;
        final int nPoints = params.nPoints;
        final int nDer = params.nDer;
        double temperatureK = params.temperature;
        long steps = params.numSteps;

        double refFrac = params.refFrac;
        double sigmaHSRef = params.sigmaHSRef;
        int[] seed = params.seed;

        boolean doHist = params.doHist;
        boolean dorefpref = params.dorefpref;
        boolean doChainRef = params.doChainRef;

        final double BDtol = params.BDtol;

        double temperature = Kelvin.UNIT.toSim(temperatureK);

        long blockSize = 1000;
        int EqSubSteps = 1000;

        double vhs = (4.0 / 3.0) * Math.PI * sigmaHSRef * sigmaHSRef * sigmaHSRef;
        final double HSBn = doChainRef ? SpecialFunctions.factorial(nPoints) / 2 * Math.pow(vhs, nPoints - 1) : Standard.BHS(nPoints, sigmaHSRef);

        System.out.println("Overlap sampling for TraPPE "+ chemForm + " at " + temperatureK + " K " + "for B"+nPoints+" and "+nDer+" derivatives");

        System.out.println("Reference diagram: B"+nPoints+" for hard spheres with diameter " + sigmaHSRef + " Angstroms");

        System.out.println("  B"+nPoints+"HS: "+HSBn);

        System.out.println(steps + " steps (" + (steps / blockSize) + " blocks of " + blockSize + ")");

        Space space = Space3D.getInstance();

        MayerFunction fRefPos = new MayerFunction() {
            public void setBox(Box box) {
            }

            public IPotential getPotential() {
                return null;
            }

            public double f(IMoleculeList pair, double r2, double beta) {
                return r2 < sigmaHSRef * sigmaHSRef ? 1 : 0;
            }
        };

        // Setting up reference cluster
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        ClusterAbstract refCluster = doChainRef ? new ClusterChainHS(nPoints, fRefPos) : new ClusterWheatleyHS(nPoints, fRef);
        refCluster.setTemperature(temperature);

        //Setting up target cluster
        PotentialGroup pTargetGroup = null;
        Species species = null;

        if(chemForm == ChemForm.N2) {

            double bondLength = 1.10; // Angstrom
            double sigmaN = 3.31; // Angstrom
            double epsilonN = Kelvin.UNIT.toSim(36.0);
            double qN = Electron.UNIT.toSim(-0.482);
            double qM = Electron.UNIT.toSim(0.964);



            Vector3D posM = new Vector3D(new double[] {0,0,0});
            Vector3D posN1 = new Vector3D(new double[] {-bondLength/2,0,0});
            Vector3D posN2 = new Vector3D(new double[] {+bondLength/2,0,0});
            Vector[] pos = new Vector[]{posM,posN1,posN2};
            IConformation conformation = new ConformationGeneric(pos);

            pTargetGroup = new PotentialGroup(2);

            P2LennardJones p2N = new P2LennardJones(space, sigmaN, epsilonN);
            P2Electrostatic p2eNN = new P2Electrostatic(space);
            p2eNN.setCharge1(qN);
            p2eNN.setCharge2(qN);
            P2Electrostatic p2eMM = new P2Electrostatic(space);
            p2eMM.setCharge1(qM);
            p2eMM.setCharge2(qM);
            P2ElectrostaticWithHardCore p2eNM = new P2ElectrostaticWithHardCore(space);
            p2eNM.setCharge1(qN);
            p2eNM.setCharge2(qM);
            p2eNM.setSigma(0.1);

            AtomType typeM = new AtomType(new ElementSimple("M", 0.0));
            AtomType typeN = new AtomType(Nitrogen.INSTANCE);

            AtomType[] atomTypes = new AtomType[]{typeM,typeN};

            species = new SpeciesSpheresHetero(space, atomTypes);

            species.setConformation(conformation);
            ((SpeciesSpheresHetero) species).setChildCount(new int[]{1,2});

            pTargetGroup.addPotential(p2N, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeN, typeN}));
//            pTargetGroup.addPotential(p2eNN, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeN, typeN}));
//            pTargetGroup.addPotential(p2eMM, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeM, typeM}));
//            pTargetGroup.addPotential(p2eNM, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeN, typeM}));

        }
        else if (chemForm == ChemForm.O2) {

            double bondLength = 1.210; // Angstrom
            double sigmaO = 3.020; // Angstrom
            double epsilonO = Kelvin.UNIT.toSim(49.0);
            double qO = Electron.UNIT.toSim(-0.113);
            double qM = Electron.UNIT.toSim(0.226);

            Vector3D posM = new Vector3D(new double[] {0,0,0});
            Vector3D posO1 = new Vector3D(new double[] {-bondLength/2,0,0});
            Vector3D posO2 = new Vector3D(new double[] {+bondLength/2,0,0});
            Vector[] pos = new Vector[]{posM,posO1,posO2};
            IConformation conformation = new ConformationGeneric(pos);

            pTargetGroup = new PotentialGroup(2);

            P2LennardJones p2O = new P2LennardJones(space, sigmaO, epsilonO);
            P2Electrostatic p2eOO = new P2Electrostatic(space);
            p2eOO.setCharge1(qO);
            p2eOO.setCharge2(qO);
            P2Electrostatic p2eMM = new P2Electrostatic(space);
            p2eMM.setCharge1(qM);
            p2eMM.setCharge2(qM);
            P2ElectrostaticWithHardCore p2eOM = new P2ElectrostaticWithHardCore(space);
            p2eOM.setCharge1(qO);
            p2eOM.setCharge2(qM);
            p2eOM.setSigma(0.1);

            AtomType typeM = new AtomType(new ElementSimple("M", 0.0));
            AtomType typeO = new AtomType(Oxygen.INSTANCE);

            AtomType[] atomTypes = new AtomType[]{typeM,typeO};

            species = new SpeciesSpheresHetero(space, atomTypes);

            species.setConformation(conformation);
            ((SpeciesSpheresHetero) species).setChildCount(new int[]{1,2});

            pTargetGroup.addPotential(p2O, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeO, typeO}));
            pTargetGroup.addPotential(p2eOO, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeO, typeO}));
            pTargetGroup.addPotential(p2eMM, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeM, typeM}));
            pTargetGroup.addPotential(p2eOM, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeO, typeM}));

        }
        else if (chemForm == ChemForm.CO2) {

            double bondLengthCO = 1.160; // Angstrom
            double sigmaC = 2.800; // Angstrom
            double epsilonC = Kelvin.UNIT.toSim(27.0);
            double qC = Electron.UNIT.toSim(0.700);
            double sigmaO = 3.050; // Angstrom
            double epsilonO = Kelvin.UNIT.toSim(79.0);
            double qO = Electron.UNIT.toSim(-0.350);

            double sigmaCO = (sigmaC+sigmaO)/2;
            double epsilonCO = Math.sqrt(epsilonC*epsilonO);

            Vector3D posC = new Vector3D(new double[] {0,0,0});
            Vector3D posO1 = new Vector3D(new double[] {-bondLengthCO,0,0});
            Vector3D posO2 = new Vector3D(new double[] {+bondLengthCO,0,0});
            Vector[] pos = new Vector[]{posC,posO1,posO2};
            IConformation conformation = new ConformationGeneric(pos);

            pTargetGroup = new PotentialGroup(2);

            P2LennardJones p2O = new P2LennardJones(space, sigmaO, epsilonO);
            P2LennardJones p2C = new P2LennardJones(space, sigmaC, epsilonC);
            P2LennardJones p2CO = new P2LennardJones(space, sigmaCO, epsilonCO);
            P2Electrostatic p2eOO = new P2Electrostatic(space);
            p2eOO.setCharge1(qO);
            p2eOO.setCharge2(qO);
            P2Electrostatic p2eCC = new P2Electrostatic(space);
            p2eCC.setCharge1(qC);
            p2eCC.setCharge2(qC);
            P2Electrostatic p2eCO = new P2Electrostatic(space);
            p2eCO.setCharge1(qO);
            p2eCO.setCharge2(qC);

            AtomType typeC = new AtomType(Carbon.INSTANCE);
            AtomType typeO = new AtomType(Oxygen.INSTANCE);

            AtomType[] atomTypes = new AtomType[]{typeC,typeO};

            species = new SpeciesSpheresHetero(space, atomTypes);

            species.setConformation(conformation);
            ((SpeciesSpheresHetero) species).setChildCount(new int[]{1,2});

            pTargetGroup.addPotential(p2O, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeO, typeO}));
            pTargetGroup.addPotential(p2C, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeC, typeC}));
            pTargetGroup.addPotential(p2CO, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeC, typeO}));
            pTargetGroup.addPotential(p2eOO, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeO, typeO}));
            pTargetGroup.addPotential(p2eCC, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeC, typeC}));
            pTargetGroup.addPotential(p2eCO, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeC, typeO}));
        }
        else if (chemForm == ChemForm.NH3) {

            double bondLengthNH = 1.012; // Angstrom
            double bondLengthNM = 0.080; // Angstrom
            double thetaHNH = Degree.UNIT.toSim(106.7);
            double thetaHNM = Degree.UNIT.toSim(67.9) ;
            double sigmaN = 3.420; // Angstrom
            double epsilonN = Kelvin.UNIT.toSim(185.0);
            double qH = Electron.UNIT.toSim(0.410);
            double qM = Electron.UNIT.toSim(-1.230);

            Vector3D posN = new Vector3D(new double[] {0,0,0});
            Vector3D posH1 = new Vector3D(new double[] {bondLengthNH,0,0});
            Vector3D posH2 = new Vector3D(new double[] {bondLengthNH*(Math.sin(thetaHNH)),-bondLengthNH*(Math.cos(thetaHNH)),0});
            Vector3D posH3 = new Vector3D(new double[] {});
            Vector3D posM = new Vector3D(new double[] {});
            Vector[] pos = new Vector[]{posN,posH1,posH2,posH3,posM};
            IConformation conformation = new ConformationGeneric(pos);

            pTargetGroup = new PotentialGroup(2);

            P2LennardJones p2N = new P2LennardJones(space, sigmaN, epsilonN);
            P2Electrostatic p2eHH = new P2Electrostatic(space);
            p2eHH.setCharge1(qH);
            p2eHH.setCharge2(qH);
            P2Electrostatic p2eMM = new P2Electrostatic(space);
            p2eMM.setCharge1(qM);
            p2eMM.setCharge2(qM);
            P2ElectrostaticWithHardCore p2eHM = new P2ElectrostaticWithHardCore(space);
            p2eHM.setCharge1(qH);
            p2eHM.setCharge2(qM);
            p2eHM.setSigma(0.1);

            AtomType typeN = new AtomType(Nitrogen.INSTANCE);
            AtomType typeH = new AtomType(Hydrogen.INSTANCE);
            AtomType typeM = new AtomType(new ElementSimple("M", 0.0));

            AtomType[] atomTypes = new AtomType[]{typeN,typeH,typeM};

            species = new SpeciesSpheresHetero(space, atomTypes);

            species.setConformation(conformation);
            ((SpeciesSpheresHetero) species).setChildCount(new int[]{1,3,1});

            pTargetGroup.addPotential(p2N, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeN, typeN}));
            pTargetGroup.addPotential(p2eHH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeH, typeH}));
            pTargetGroup.addPotential(p2eMM, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeM, typeM}));
            pTargetGroup.addPotential(p2eHM, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeH, typeM}));
        }

        MayerGeneral fTarget = new MayerGeneral(pTargetGroup);

        final ClusterWheatleySoftDerivatives targetCluster = new ClusterWheatleySoftDerivatives(nPoints, fTarget, BDtol, nDer);
        targetCluster.setTemperature(temperature);

        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, species, nPoints, temperature,refCluster,targetCluster);
        if(seed!=null)sim.setRandom(new RandomMersenneTwister(seed));
        System.out.println("random seeds: "+ Arrays.toString(seed==null?sim.getRandomSeeds():seed));

        ClusterMultiToSingle[] primes = new ClusterMultiToSingle[nDer];
        for(int m=0;m<primes.length;m++){
            primes[m]= new ClusterMultiToSingle(targetCluster, m+1);
        }
        sim.setExtraTargetClusters(primes);

        sim.init();

        if (doChainRef) {
            sim.integrators[0].getMoveManager().removeMCMove(sim.mcMoveTranslate[0]);
            MCMoveClusterMoleculeHSChain mcMoveHSC = new MCMoveClusterMoleculeHSChain(sim.getRandom(), space, sigmaHSRef);
            sim.integrators[0].getMoveManager().addMCMove(mcMoveHSC);
            sim.accumulators[0].setBlockSize(1);
        }

        sim.integratorOS.setNumSubSteps(EqSubSteps);

        sim.integratorOS.setAggressiveAdjustStepFraction(true);

        if (false) {
            sim.box[0].getBoundary().setBoxSize(space.makeVector(new double[]{10,10,10}));
            sim.box[1].getBoundary().setBoxSize(space.makeVector(new double[]{10,10,10}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box[0]);
            DisplayBox displayBox1 = simGraphic.getDisplayBox(sim.box[1]);
//            displayBox0.setPixelUnit(new Pixel(300.0/size));
//            displayBox1.setPixelUnit(new Pixel(300.0/size));
            displayBox0.setShowBoundary(false);
            displayBox1.setShowBoundary(false);
            ((DisplayBoxCanvasG3DSys)displayBox0.canvas).setBackgroundColor(Color.WHITE);
            ((DisplayBoxCanvasG3DSys)displayBox1.canvas).setBackgroundColor(Color.WHITE);


//            ColorSchemeRandomByMolecule colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box[0], sim.getRandom());
//            displayBox0.setColorScheme(colorScheme);
//            colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box[1], sim.getRandom());
//            displayBox1.setColorScheme(colorScheme);
            simGraphic.makeAndDisplayFrame();
            DiameterHashByType dhbt = (DiameterHashByType)simGraphic.getDisplayBox(sim.box[0]).getDiameterHash();
            dhbt.setDiameter(species.getAtomType(0), 0.5);
            dhbt.setDiameter(species.getAtomType(1), 0.5);
            simGraphic.getDisplayBox(sim.box[1]).setDiameterHash(dhbt);

            ColorSchemeByType csbt = (ColorSchemeByType)simGraphic.getDisplayBox(sim.box[0]).getColorScheme();
            csbt.setColor(species.getAtomType(0), Color.BLACK);
            simGraphic.getDisplayBox(sim.box[1]).setColorScheme(csbt);

            sim.integratorOS.setNumSubSteps(1000);
            sim.setAccumulatorBlockSize(1000);

            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
            sim.getController().removeAction(sim.ai);
            sim.getController().addAction(new IAction() {
                public void actionPerformed() {
                    sim.initRefPref(null, 10);
                    sim.equilibrate(null, 20);
                    sim.ai.setMaxSteps(Long.MAX_VALUE);
                }
            });
            sim.getController().addAction(sim.ai);
            if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
                throw new RuntimeException("Oops");
            }
            return;
        }

        long t1 = System.currentTimeMillis();

        if (refFrac >= 0) {
            sim.integratorOS.setRefStepFraction(refFrac);
            sim.integratorOS.setAdjustStepFraction(false);
        }

        System.out.println();
        String refFileName = null;

        if (isCommandline) {
            // if running interactively, don't use the file
            String tempString = ""+temperatureK;
            if (temperatureK == (int)temperatureK) {
                // temperature is an integer, use "200" instead of "200.0"
                tempString = ""+(int)temperatureK;
            }
            refFileName = "refpref_"+"_"+nPoints+"_"+tempString+"K";
        }


        final HistogramSimple targHist = new HistogramSimple(200, new DoubleRange(-1, 4));
        IntegratorListener histListenerTarget = new IntegratorListener() {
            public void integratorStepStarted(IntegratorEvent e) {}

            public void integratorStepFinished(IntegratorEvent e) {
                CoordinatePairSet cPairs = sim.box[1].getCPairSet();
                for (int i=0; i<nPoints; i++) {
                    for (int j=i+1; j<nPoints; j++) {
                        double r2 = cPairs.getr2(i, j);
                        double r = Math.sqrt(r2);
                        if (r > 1) {
                            r = Math.log(r);
                        }
                        else {
                            r -= 1;
                        }
                        targHist.addValue(r);
                    }
                }

            }

            public void integratorInitialized(IntegratorEvent e) {}
        };

        if (doHist) {
            System.out.println("collecting histograms");
            // only collect the histogram if we're forcing it to run the reference system
            sim.integrators[1].getEventManager().addListener(histListenerTarget);
        }

        sim.initRefPref(refFileName, (steps / EqSubSteps) / 20);
        sim.equilibrate(refFileName, (steps / EqSubSteps) / 10);

        System.out.println("equilibration finished");

        if(dorefpref){
            long t2 = System.currentTimeMillis();
            System.out.println("time: "+(t2-t1)/1000.0);
            return;
        }

        IntegratorListener progressReport = new IntegratorListener() {

            public void integratorStepStarted(IntegratorEvent e) {}

            public void integratorStepFinished(IntegratorEvent e) {
                if (sim.integratorOS.getStepCount() % 100 != 0) return;
                System.out.print(sim.integratorOS.getStepCount()+" steps: ");
                double[] ratioAndError = sim.dvo.getAverageAndError();
                System.out.println("abs average: "+ratioAndError[0]*HSBn+", error: "+ratioAndError[1]*HSBn);
            }

            public void integratorInitialized(IntegratorEvent e) {}

        };
        if (false) {
            sim.integratorOS.getEventManager().addListener(progressReport);
        }

        sim.integratorOS.setNumSubSteps((int) blockSize);
        sim.setAccumulatorBlockSize(blockSize);

        if (doChainRef) sim.accumulators[0].setBlockSize(1);
        sim.ai.setMaxSteps(steps / blockSize);
        for (int i=0; i<2; i++) {
            if (i > 0 || !doChainRef) System.out.println("MC Move step sizes " + sim.mcMoveTranslate[i].getStepSize());
        }

        sim.getController().actionPerformed();

        if (doHist) {
            double[] xValues = targHist.xValues();
            double[] h = targHist.getHistogram();
            for (int i=0; i<xValues.length; i++) {
                if (!Double.isNaN(h[i])) {
                    double r = xValues[i];
                    double y = h[i];
                    if (r < 0) r += 1;
                    else {
                        r = Math.exp(r);
                        y /= r;
                    }
                    System.out.println(r+" "+y);
                }
            }
        }

        System.out.println("final reference step fraction "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step fraction "+sim.integratorOS.getRefStepFraction());

        String[] extraNames = new String[nDer];
        for (int i = 1; i <= nDer; i++) {
            extraNames[i - 1] = "derivative " + i;
        }
        sim.printResults(HSBn, extraNames);

        long t2 = System.currentTimeMillis();
        System.out.println("time: "+(t2-t1)/1000.0);

	}

    enum ChemForm {
        N2, O2, CO2, NH3
    }

    /**
     * Inner class for parameters
     */
    public static class VirialParam extends ParameterBase {
        // don't change these
        public ChemForm chemForm = ChemForm.N2;
        public int nPoints = 2;
        public int nDer = 3;
        public double temperature = 400;
        public long numSteps = 1000000;
        public double refFrac = -1;
        public double sigmaHSRef = 5;
        public boolean doHist = false;
        public int[] seed = null;
        public boolean dorefpref = false;
        public boolean doChainRef = true;
        public double BDtol = 1e-12;

    }
}
