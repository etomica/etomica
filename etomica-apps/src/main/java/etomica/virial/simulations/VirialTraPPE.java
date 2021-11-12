/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.*;
import etomica.graphics.ColorSchemeRandomByMolecule;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.SimulationGraphic;
import etomica.math.SpecialFunctions;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculePositionCOM;
import etomica.potential.IPotential;
import etomica.potential.P2PotentialGroupBuilder;
import etomica.potential.PotentialMoleculePair;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesManager;
import etomica.units.Degree;
import etomica.units.Electron;
import etomica.units.Kelvin;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;
import etomica.virial.MayerFunction;
import etomica.virial.MayerGeneral;
import etomica.virial.MayerHardSphere;
import etomica.virial.cluster.*;
import etomica.virial.mcmove.MCMoveClusterMoleculeHSChain;
import etomica.virial.wheatley.ClusterWheatleyHS;
import etomica.virial.wheatley.ClusterWheatleySoftDerivatives;
import etomica.virial.wheatley.ClusterWheatleySoftDerivativesMix;
import etomica.virial.wheatley.ClusterWheatleySoftDerivativesMixBD;

import java.awt.*;
import java.util.Arrays;

/**
 * Compute pure, binary, ternary and quaternary mixture virial coefficients using overlap sampling simulations
 * for some molecules using the TraPPE force fields.
 *
 */
public class VirialTraPPE {


    public static void main(String[] args) {
        VirialParam params = new VirialParam();
        boolean isCommandline = args.length > 0;
        if (isCommandline) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            // Customize Interactive Parameters Here
            params.chemForm = new ChemForm[]{ChemForm.NH3};
            params.nPoints = 6;
            params.nTypes = new int[]{6};
            params.nDer = 3;
            params.temperature = 530;
            params.numSteps = 1000000;

            params.refFrac = -1;
            params.sigmaHSRef = 5;
            params.seed = null;

            params.dorefpref = false;
            params.doChainRef = true;

            params.BDtol = 1e-12;
        }

        // Import Params
        final ChemForm[] chemForm = params.chemForm;
        final int nPoints = params.nPoints;
        final int[] nTypes = params.nTypes;
        final int nDer = params.nDer;
        final double temperatureK = params.temperature;
        final long steps = params.numSteps;

        double refFrac = params.refFrac;
        double sigmaHSRef = params.sigmaHSRef;
        int[] seed = params.seed;

        boolean dorefpref = params.dorefpref;
        boolean doChainRef = params.doChainRef;

        double BDtol = params.BDtol;


        // Set Number of Blocks
        final long numBlocks = 1000;

        // Set Big Decimal Acceptance Fraction
        final double BDAccFrac = 0.1;

        // Check Params
        if( chemForm.length != nTypes.length ) throw new RuntimeException("chemFrom and nTypes lengths are unequal!");

        if( chemForm.length > 1 && Arrays.stream(nTypes).sum() != nPoints ) throw new RuntimeException("nPoints and nTypes do not match!");

        boolean isMixture = ( nTypes.length > 1 ) ;

        // Check if Pure or Mixture
        if(isMixture){
            for(int i=0; i<nTypes.length; i++){
                if(nTypes[i]==nPoints) isMixture=false;
            }
        }

        // Set Simulation Temperature
        double temperature = Kelvin.UNIT.toSim(temperatureK);

        // Evaluate Hard Sphere Coefficient
        double vhs = (4.0 / 3.0) * Math.PI * sigmaHSRef * sigmaHSRef * sigmaHSRef;
        final double HSBn = doChainRef ? SpecialFunctions.factorial(nPoints) / 2 * Math.pow(vhs, nPoints - 1) : Standard.BHS(nPoints, sigmaHSRef);

        // Print Pretext
        if(!isMixture) {
            ChemForm chemFormPure = chemForm[0];
            if(nTypes.length>1) {
                for(int i=0; i<nTypes.length; i++){
                    if(nTypes[i]==nPoints) chemFormPure=chemForm[i];
                }
            }

            System.out.println("Overlap sampling for TraPPE " + chemFormPure + " at " + temperatureK + " K " + "for B" + nPoints + " and " + nDer + " derivatives");
        }
        else{
            String nTstr="{";
            for(int i=0; i<nTypes.length; i++){
                if(nTypes[i]!=0) nTstr += ((nTstr=="{") ? "":",")+nTypes[i];
            }
            nTstr+="}";

            String CFstr="";
            for(int i=0; i<chemForm.length; i++){
                if(nTypes[i]!=0) CFstr += chemForm[i]+" ";
            }

            System.out.println("Overlap sampling for TraPPE " + CFstr + " " +nTstr + " Mixture at " + temperatureK + " K " + "for B" + nPoints + " and " + nDer + " derivatives");
        }

        System.out.println("Reference diagram: B"+nPoints+" for hard spheres with diameter " + sigmaHSRef + " Angstroms");

        System.out.println("  B"+nPoints+"HS: "+HSBn);

        // Set up Space
        Space space = Space3D.getInstance();

        // Setting up Reference Cluster Mayer Function
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

        // Setting up Reference Cluster
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        ClusterAbstract refCluster = doChainRef ? new ClusterChainHS(nPoints, fRefPos) : new ClusterWheatleyHS(nPoints, fRef);
        refCluster.setTemperature(temperature);

        // Setting up Target Cluster Mayer Function
        ISpecies[] species = null;
        ClusterAbstractMultivalue targetCluster = null;
        ClusterWheatleySoftDerivativesMixBD targetClusterBD = null;
        SpeciesManager.Builder sb = SpeciesManager.builder();

        boolean anyPolar = false;
        MayerFunction[][] fAll = new MayerFunction[nTypes.length][nTypes.length];
        species = new ISpecies[chemForm.length];

        TraPPEParams[] TPList = new TraPPEParams[chemForm.length];

        for(int i=0; i<TPList.length; i++){
            TPList[i] = new TraPPEParams(space, chemForm[i]);
        }

        for(int i=0; i<chemForm.length; i++) {
            sb.addSpecies(TPList[i].species);
        }
        SpeciesManager sm = sb.build();

        for(int i=0; i<chemForm.length; i++) {
            TraPPEParams TPi = TPList[i];
            TPi.buildPotentials(sm);
            PotentialMoleculePair PGii = TPi.potentialGroup;

            P2PotentialGroupBuilder.ModelParams MPi = new P2PotentialGroupBuilder.ModelParams(TPi.atomTypes,TPi.sigma,TPi.epsilon,TPi.charge);
            fAll[i][i] = new MayerGeneral(PGii);

            anyPolar=(anyPolar||TPi.polar);

            for(int j=i+1; j<chemForm.length; j++){

                TraPPEParams TPj = TPList[j];

                P2PotentialGroupBuilder.ModelParams MPj = new P2PotentialGroupBuilder.ModelParams(TPj.atomTypes,TPj.sigma,TPj.epsilon,TPj.charge);

                PotentialMoleculePair PGij = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space,sm, MPi,MPj);

                fAll[i][j] = fAll[j][i] = new MayerGeneral(PGij);

            }
        }

        // Setting up Target Cluster
        targetCluster = new ClusterWheatleySoftDerivativesMix(nPoints, nTypes,fAll, BDtol, nDer);
        targetCluster.setTemperature(temperature);


        // Setting BlockSize
        long blockSize = steps/numBlocks;
        int EqSubSteps = 1000;

        System.out.println(steps + " steps (" + numBlocks + " blocks of " + blockSize + ")");
        System.out.println("BD_Tol: " + BDtol + " BDAccFrac: " + BDAccFrac);

        // Setting up Flipping
        if(anyPolar && nPoints==2) {
            System.out.println("Performing Flipping");
            ((ClusterWheatleySoftDerivativesMix) targetCluster).setTolerance(0);
            final int precision = -3*(int)Math.log10(BDtol);
            targetClusterBD = new ClusterWheatleySoftDerivativesMixBD(nPoints,nTypes,fAll,precision,nDer);
            targetClusterBD.setTemperature(temperature);
            ((ClusterWheatleySoftDerivativesMix) targetCluster).setDoCaching(false);
            targetClusterBD.setDoCaching(false);
            targetCluster = new ClusterCoupledFlippedMultivalue(targetCluster, targetClusterBD, space, 20, nDer, BDtol);
        }

        // Setting up Simulation
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, species, nTypes, temperature,refCluster,targetCluster);
        sim.setDoFasterer(true);
        if(seed!=null)sim.setRandom(new RandomMersenneTwister(seed));
        System.out.println("random seeds: "+ Arrays.toString(seed==null?sim.getRandomSeeds():seed));
        if(targetCluster instanceof ClusterCoupledFlippedMultivalue) {
            ((ClusterCoupledFlippedMultivalue) targetCluster).setBDAccFrac(BDAccFrac,sim.getRandom());
        }
        else {
            ((ClusterWheatleySoftDerivativesMix) targetCluster).setBDAccFrac(BDAccFrac, sim.getRandom());
            ((ClusterWheatleySoftDerivativesMix) targetCluster).setNumBDCheckBins(8);
        }

        // Adding derivative clusters to simulation
        ClusterMultiToSingle[] primes = new ClusterMultiToSingle[nDer];
        for(int m=0;m<primes.length;m++){
            primes[m] = new ClusterMultiToSingle(targetCluster, m + 1);
        }
        sim.setExtraTargetClusters(primes);

        // Initialize Simulation
        sim.init();

        // Set Position Definitions
        sim.box[0].setPositionDefinition(new MoleculePositionCOM(space));
        sim.box[1].setPositionDefinition(new MoleculePositionCOM(space));

        // Setting Chain Ref Moves
        if (doChainRef) {
            sim.integratorsFasterer[0].getMoveManager().removeMCMove(sim.mcMoveTranslate[0]);
            MCMoveClusterMoleculeHSChain mcMoveHSC = new MCMoveClusterMoleculeHSChain(sim.getRandom(), sim.box[0], sigmaHSRef);
            sim.integratorsFasterer[0].getMoveManager().addMCMove(mcMoveHSC);
            sim.accumulators[0].setBlockSize(1);
        }

        // Run with Graphics
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


            ColorSchemeRandomByMolecule colorScheme = new ColorSchemeRandomByMolecule(sim.getSpeciesManager(), sim.box[0], sim.getRandom());
            displayBox0.setColorScheme(colorScheme);
            colorScheme = new ColorSchemeRandomByMolecule(sim.getSpeciesManager(), sim.box[1], sim.getRandom());
            displayBox1.setColorScheme(colorScheme);
            simGraphic.makeAndDisplayFrame();

            sim.integratorOS.setNumSubSteps(1000);
            sim.setAccumulatorBlockSize(1000);

            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
            sim.initRefPref(null, 10, false);
            sim.equilibrate(null, 20, false);
            sim.getController().addActivity(new ActivityIntegrate(sim.integratorOS));
            return;
        }

        // Setting up Equilibration
        sim.integratorOS.setNumSubSteps(EqSubSteps);
        sim.integratorOS.setAggressiveAdjustStepFraction(true);

        // Start timing
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

        // Equilibrate
        sim.initRefPref(refFileName, (steps / EqSubSteps) / 20);
        sim.equilibrate(refFileName, (steps / EqSubSteps) / 10);

        System.out.println("equilibration finished");

        if(dorefpref){
            long t2 = System.currentTimeMillis();
            System.out.println("time: "+(t2-t1)/1000.0);
            return;
        }

        // Setting up Production Run
        sim.integratorOS.setNumSubSteps((int) blockSize);
        sim.setAccumulatorBlockSize(blockSize);

        if (doChainRef) sim.accumulators[0].setBlockSize(1);
        for (int i = 0; i < 2; i++) {
            if (i > 0 || !doChainRef) System.out.println("MC Move step sizes " + sim.mcMoveTranslate[i].getStepSize());
        }

        // Production Run
        ActivityIntegrate ai = new ActivityIntegrate(sim.integratorOS, numBlocks);
        sim.getController().runActivityBlocking(ai);

        System.out.println();

        // Print Simulation Output
        System.out.println("final reference step fraction " + sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step fraction " + sim.integratorOS.getRefStepFraction());

        String[] extraNames = new String[nDer];
        for (int i = 1; i <= nDer; i++) {
            extraNames[i - 1] = "derivative " + i;
        }
        sim.printResults(HSBn, extraNames);

        // Print BD and Flip Stats
        if (targetCluster instanceof ClusterWheatleySoftDerivatives) {
            System.out.println("SoftBDcount: " + ((ClusterWheatleySoftDerivatives)targetCluster).getSoftBDcount() + " SoftBDfrac: " + ((ClusterWheatleySoftDerivatives)targetCluster).getSoftBDfrac() + " Softcount: " + ((ClusterWheatleySoftDerivatives)targetCluster).getSoftcount());
        }
        else if (targetCluster instanceof ClusterCoupledFlippedMultivalue) {
            ClusterCoupledFlippedMultivalue foo = (ClusterCoupledFlippedMultivalue)targetCluster;
            System.out.println("BDcount: " + foo.getBDcount() + " BDfrac: " + foo.getBDfrac() + " totBDcount: " + foo.getBDtotcount());
            System.out.println("FlipCount: " + foo.getflipcount() + " Flipfrac: " + foo.getflipfrac() + " FlipTotcount: " + foo.gettotcount());
        }

        long t2 = System.currentTimeMillis();
        System.out.println("time: " + (t2 - t1) / 1000.0);

    }

    enum ChemForm {
        N2, O2, CO2, NH3
    }

    /**
     * Inner class for parameters
     */
    public static class VirialParam extends ParameterBase {
        // don't change these
        public ChemForm[] chemForm = {ChemForm.N2};
        public int nPoints = 2;
        public int[] nTypes = {2};
        public int nDer = 3;
        public double temperature = 400;
        public long numSteps = 1000000;

        public double refFrac = -1;
        public double sigmaHSRef = 5;
        public int[] seed = null;

        public boolean dorefpref = false;
        public boolean doChainRef = true;

        public double BDtol = 1e-12;

    }

    /**
     * Inner class for TraPPE Parameters
     */
    public static class TraPPEParams{

        public final AtomType[] atomTypes;
        public final double[] sigma;
        public final double[] epsilon;
        public final double[] charge;
        public final ISpecies species;
        public PotentialMoleculePair potentialGroup;
        protected static Element elementM = new ElementSimple("M", 0.0);
        protected boolean polar;
        protected final ChemForm chemForm;
        protected final Space space;
        //Set up computing the boolean. It is hard coded for now.

        public TraPPEParams(Space space, ChemForm chemForm){
            this.chemForm = chemForm;
            this.space = space;

            if(chemForm == ChemForm.N2) {

                //Atoms in Compound
                AtomType typeM = new AtomType(elementM);
                AtomType typeN = new AtomType(Nitrogen.INSTANCE);

                atomTypes = new AtomType[]{typeM,typeN};

                //TraPPE Parameters
                double bondLength = 1.10; // Angstrom
                double sigmaN = 3.31; // Angstrom
                double epsilonN = Kelvin.UNIT.toSim(36.0);
                double qN = Electron.UNIT.toSim(-0.482);
                double sigmaM = 0.0; // Angstrom
                double epsilonM = Kelvin.UNIT.toSim(0.0);
                double qM = Electron.UNIT.toSim(0.964);

                //Construct Arrays
                sigma = new double[]{sigmaM, sigmaN};
                epsilon = new double[]{epsilonM, epsilonN};
                charge = new double[]{qM, qN};

                //Get Coordinates
                Vector3D posM = new Vector3D(new double[]{0, 0, 0});
                Vector3D posN1 = new Vector3D(new double[]{-bondLength / 2, 0, 0});
                Vector3D posN2 = new Vector3D(new double[]{+bondLength / 2, 0, 0});

                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeM, posM, "M")
                        .addAtom(typeN, posN1, "N1")
                        .addAtom(typeN, posN2, "N2")
                        .build();
            }

            else if (chemForm == ChemForm.O2) {

                //Atoms in Compound
                AtomType typeM = new AtomType(elementM);
                AtomType typeO = new AtomType(Oxygen.INSTANCE);

                atomTypes = new AtomType[]{typeM,typeO};

                //TraPPE Parameters
                double bondLength = 1.210; // Angstrom
                double sigmaO = 3.020; // Angstrom
                double epsilonO = Kelvin.UNIT.toSim(49.0);
                double qO = Electron.UNIT.toSim(-0.113);
                double sigmaM = 0.0; // Angstrom
                double epsilonM = Kelvin.UNIT.toSim(0.0);
                double qM = Electron.UNIT.toSim(0.226);

                //Construct Arrays
                sigma = new double[]{sigmaM, sigmaO};
                epsilon = new double[]{epsilonM, epsilonO};
                charge = new double[]{qM, qO};

                //Get Coordinates
                Vector3D posM = new Vector3D(new double[]{0, 0, 0});
                Vector3D posO1 = new Vector3D(new double[]{-bondLength / 2, 0, 0});
                Vector3D posO2 = new Vector3D(new double[]{+bondLength / 2, 0, 0});

                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeM, posM, "M")
                        .addAtom(typeO, posO1, "O1")
                        .addAtom(typeO, posO2, "O2")
                        .build();
            }

            else if (chemForm == ChemForm.CO2) {

                //Atoms in Compound
                AtomType typeC = new AtomType(Carbon.INSTANCE);
                AtomType typeO = new AtomType(Oxygen.INSTANCE);

                atomTypes = new AtomType[]{typeC,typeO};

                //TraPPE Parameters
                double bondLengthCO = 1.160; // Angstrom
                double sigmaC = 2.800; // Angstrom
                double epsilonC = Kelvin.UNIT.toSim(27.0);
                double qC = Electron.UNIT.toSim(0.700);
                double sigmaO = 3.050; // Angstrom
                double epsilonO = Kelvin.UNIT.toSim(79.0);
                double qO = Electron.UNIT.toSim(-0.350);

                //Construct Arrays
                sigma = new double[]{sigmaC, sigmaO};
                epsilon = new double[]{epsilonC, epsilonO};
                charge = new double[]{qC, qO};

                //Get Coordinates
                Vector3D posC = new Vector3D(new double[]{0, 0, 0});
                Vector3D posO1 = new Vector3D(new double[]{-bondLengthCO, 0, 0});
                Vector3D posO2 = new Vector3D(new double[]{+bondLengthCO, 0, 0});

                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeC, posC, "C")
                        .addAtom(typeO, posO1, "O1")
                        .addAtom(typeO, posO2, "O2")
                        .build();
            }
            else if (chemForm == ChemForm.NH3) {

                //Atom in Compound
                AtomType typeN = new AtomType(Nitrogen.INSTANCE);
                AtomType typeH = new AtomType(Hydrogen.INSTANCE);
                AtomType typeM = new AtomType(elementM);

                atomTypes = new AtomType[]{typeN,typeH,typeM};

                polar = true;

                //TraPPE Parameters
                double bondLengthNH = 1.012; // Angstrom
                double bondLengthNM = 0.080; // Angstrom
                double thetaHNM = Degree.UNIT.toSim(67.9) ;
                double thetaHNH = Degree.UNIT.toSim(106.7);
                double thetaHNHxy = Degree.UNIT.toSim(60);
                double sigmaN = 3.420; // Angstrom
                double epsilonN = Kelvin.UNIT.toSim(185.0);
                double qN = Electron.UNIT.toSim(0.0);
                double sigmaH = 0.0; // Angstrom
                double epsilonH = Kelvin.UNIT.toSim(0.0);
                double qH = Electron.UNIT.toSim(0.410);
                double sigmaM = 0.0; // Angstrom
                double epsilonM = Kelvin.UNIT.toSim(0.0);
                double qM = Electron.UNIT.toSim(-1.230);

                //Construct Arrays
                sigma = new double[] {sigmaN,sigmaH,sigmaM};
                epsilon = new double[] {epsilonN,epsilonH,epsilonM};
                charge = new double[]{qN, qH, qM};

                //Get Coordinates
                Vector3D posN = new Vector3D(new double[]{0, 0, 0});
                Vector3D posH1 = new Vector3D(new double[]{bondLengthNH * Math.sin(thetaHNM), 0, -bondLengthNH * Math.cos(thetaHNM)});
                Vector3D posH2 = new Vector3D(new double[]{-bondLengthNH * Math.sin(thetaHNM) * Math.cos(thetaHNHxy), bondLengthNH * Math.sin(thetaHNM) * Math.sin(thetaHNHxy), -bondLengthNH * Math.cos(thetaHNM)});
                Vector3D posH3 = new Vector3D(new double[]{-bondLengthNH * Math.sin(thetaHNM) * Math.cos(thetaHNHxy), -bondLengthNH * Math.sin(thetaHNM) * Math.sin(thetaHNHxy), -bondLengthNH * Math.cos(thetaHNM)});
                Vector3D posM = new Vector3D(new double[]{0, 0, -bondLengthNM});

                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeN, posN, "N")
                        .addAtom(typeH, posH1, "H1")
                        .addAtom(typeH, posH2, "H2")
                        .addAtom(typeH, posH3, "H3")
                        .addAtom(typeM, posM, "M")
                        .build();
            }
            else {
                throw new RuntimeException("unrecognized chem form");
            }

        }

        public void buildPotentials(SpeciesManager sm) {
            if(chemForm == ChemForm.N2) {
                P2PotentialGroupBuilder.ModelParams modelParams = new P2PotentialGroupBuilder.ModelParams(atomTypes, sigma, epsilon, charge);
                potentialGroup = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space, sm, modelParams, null);
            }
            else if (chemForm == ChemForm.O2) {
                P2PotentialGroupBuilder.ModelParams modelParams = new P2PotentialGroupBuilder.ModelParams(atomTypes, sigma, epsilon, charge);
                potentialGroup = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space, sm, modelParams, null);
            }
            else if (chemForm == ChemForm.CO2) {
                P2PotentialGroupBuilder.ModelParams modelParams = new P2PotentialGroupBuilder.ModelParams(atomTypes, sigma, epsilon, charge);
                potentialGroup = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space, sm, modelParams, null);
            }
            else if (chemForm == ChemForm.NH3) {
                P2PotentialGroupBuilder.ModelParams modelParams = new P2PotentialGroupBuilder.ModelParams(atomTypes, sigma, epsilon, charge);
                potentialGroup = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space, sm, modelParams, null);
            }

        }

    }

}
