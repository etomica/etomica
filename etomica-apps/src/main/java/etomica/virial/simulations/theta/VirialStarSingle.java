/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations.theta;

import etomica.action.XYZWriter;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import java.io.FileWriter;

import etomica.atom.DiameterHashByType;
import etomica.atom.IAtom;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountSteps;
import etomica.data.meter.MeterRadiusGyration;
import etomica.graphics.*;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.potential.IPotential2;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.PotentialMoleculePair;
import etomica.potential.compute.*;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesManager;
import etomica.starpolymer.ConformationStarPolymerGraft;
import etomica.units.Pixel;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.collections.IntArrayList;
import etomica.util.random.RandomMersenneTwister;
import etomica.virial.mcmove.MCMoveClusterAngle;
import etomica.virial.mcmove.MCMoveClusterStretch;
import etomica.virial.simulations.theta.MCMoveClusterShuffle;

import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;



public class VirialStarSingle {

    public static void main(String[] args) {
        VirialStarParams params = new VirialStarParams();

       if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.armLength = 32;
            params.numArms = 20;
            params.temperature = 1.5;
            params.numSteps = 10000;
            params.coreSigma = 2.0;  //new
            params.epsCore = 1.0;// new
            params.ideal = false;
        }
        int numArms = params.numArms;
        int armLength = params.armLength;
        double temperature = params.temperature;
        long steps = params.numSteps;
        boolean ideal = params.ideal;
        if (armLength < 2){
            throw new RuntimeException("arm length must be at least 2");
        } //new
        Space space = Space3D.getInstance();

        ConformationStarPolymerGraft conf = new ConformationStarPolymerGraft(Space3D.getInstance(), numArms, armLength);
        conf.setSigma0(params.coreSigma); //changed
        AtomType typeCore = AtomType.simple("Core");
        AtomType typeMono = AtomType.simple("Mono"); //added
        ISpecies species = new SpeciesBuilder(Space3D.getInstance())
                .addCount(typeCore, 1)//new
                .addCount(typeMono, numArms*armLength)
                .withConformation(conf)
                .build();
        SpeciesManager sm = new SpeciesManager.Builder().addSpecies(species).build();

        PotentialMoleculePair pTarget = new PotentialMoleculePair(space, sm);
        System.out.println(numArms+" arms of length "+armLength+" at T = "+temperature);
        double sigC = params.coreSigma;   // new
        double epsC = params.epsCore;     // new
        double sigM = 1;   // mono σ
        double epsM = 1;     // mono ε

        double sigCM = 0.5*(sigC + sigM);           // Lorentz
        double epsCM = Math.sqrt(epsC * epsM);      // Berthelot

        IPotential2 p2CoreCore = new P2LennardJones(epsC,  sigC);//new
        IPotential2 p2MonoMono = new P2LennardJones(epsM,  sigM);//new
        IPotential2 p2CoreMono = new P2LennardJones(epsCM, sigCM);//new

        PotentialMasterBonding.FullBondingInfo bondingInfo = new PotentialMasterBonding.FullBondingInfo(sm) {
            @Override
            public boolean skipBondedPair(boolean isPureAtoms, IAtom iAtom, IAtom jAtom) {
                return false;
            }
        };

        // we need to do this to convince the system that the molecules are not rigid
        // if bondingInfo thinks molecules are rigid then intramolecular LJ will not be computed
        IPotential2 pBonding = new IPotential2() {
            @Override
            public double getRange() { return 2; }
            @Override
            public void u012add(double r2, double[] u012) { }
        };
        List<int[]> pairs = new ArrayList<>();
        for (int i=0; i<numArms; i++) {
            pairs.add(new int[]{0,1+armLength*i});
            for (int j = 0; j < armLength-1; j++) {
                pairs.add(new int[]{1+armLength*i+j, 1+armLength*i+j+1});
            }
        }
        bondingInfo.setBondingPotentialPair(species, pBonding, pairs);


        Simulation sim = new Simulation(Space3D.getInstance(), sm);
       //debug //sim.setRandom(new RandomMersenneTwister(1));
        sim.makeBox(new BoundaryRectangularNonperiodic(space));
        sim.box().setNMolecules(species, 1);

        System.out.println("random seeds: "+ Arrays.toString(sim.getRandomSeeds()));

        //pTarget.setAtomPotential(type, type, p2);
        //new
        pTarget.setAtomPotential(typeCore, typeCore, p2CoreCore);
        pTarget.setAtomPotential(typeMono, typeMono, p2MonoMono);
        pTarget.setAtomPotential(typeCore, typeMono, p2CoreMono);
        pTarget.setAtomPotential(typeMono, typeCore, p2CoreMono);

        PotentialMasterBonding pmBonding = new PotentialMasterBonding(sm, sim.box(), bondingInfo);
        NeighborManager nbrManager = new NeighborManagerIntra(sim.box(), bondingInfo);
        if (ideal) {
            nbrManager = new NeighborManagerIntra(sim.box(), bondingInfo) {
                public NeighborIterator makeNeighborIterator() {
                    return new NeighborIterator() {
                        @Override
                        public void iterUpNeighbors(int i, NeighborConsumer consumer) {}

                        @Override
                        public void iterDownNeighbors(int i, NeighborConsumer consumer) {}

                        @Override
                        public void iterAllNeighbors(int i, NeighborConsumer consumer) {}
                    };
                }

            };
        }
        PotentialComputePair pcPair = new PotentialComputePair(sm, sim.box(), nbrManager, pTarget.getAtomPotentials());
        PotentialComputeAggregate pc = new PotentialComputeAggregate(pmBonding, pcPair);

        IntegratorMC integrator = new IntegratorMC(pc, sim.getRandom(), temperature, sim.box());

        System.out.println(steps+" steps");

        //MCMoveClusterStretch stretchMoves = null;


        IntArrayList[] bonding = new IntArrayList[1+numArms*armLength];
        bonding[0] = new IntArrayList(numArms);

        for (int i=0; i<numArms; i++) {
            bonding[0].add(1+i*armLength);
            if (armLength > 1) {
                bonding[1 + i * armLength] = new IntArrayList(new int[]{0, 1 + i * armLength + 1});
            }
            else {
                bonding[1 + i * armLength] = new IntArrayList(new int[]{0});
            }
            for (int j = 1; j < armLength-1; j++) {
                bonding[1+armLength*i+j] = new IntArrayList(new int []{1+armLength*i+j-1, 1+armLength*i+j+1});
            }
            if (armLength > 1) {
                bonding[1 + (i+1) * armLength - 1] = new IntArrayList(new int[]{1 + (i+1) * armLength - 2});
            }
        }

        //MCMoveClusterAngle angleMove = new MCMoveClusterAngle(pc, space, bonding, sim.getRandom(), 1);
        //new
        MCMoveClusterAngle angleMove = new MCMoveClusterAngle(pc, space, bonding, sim.getRandom(), 1,
                new MCMoveClusterAngle.AtomChooserStarfl(sim.getRandom(), numArms, armLength)); angleMove.setBox(sim.box());
        integrator.getMoveManager().addMCMove(angleMove);

        MCMoveClusterShuffle shuffleMove = null;
        if (armLength > 6 ) {
            shuffleMove = new MCMoveClusterShuffle(pc, space, sim.getRandom());
            shuffleMove.setExcludedAtom(0);
            shuffleMove.setBox(sim.box());
            shuffleMove.setBonding(bonding);
            shuffleMove.setStepSizeMax(armLength - 2);
            integrator.getMoveManager().addMCMove(shuffleMove);
            ((MCMoveStepTracker) shuffleMove.getTracker()).setAcceptanceTarget(0.3);

        }


        if (false) {
    double size = (2*armLength + 5) * 1.5;
    sim.box().getBoundary().setBoxSize(Vector.of(new double[]{size, size, size}));
    sim.getController().addActivity(new ActivityIntegrate(integrator), Long.MAX_VALUE, 10);
    SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "Star Single", 1);
    DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box());
    displayBox0.setPixelUnit(new Pixel(300.0 / size));
    displayBox0.setShowBoundary(false);
    ((DisplayBoxCanvasG3DSys) displayBox0.canvas).setBackgroundColor(Color.WHITE);
    DiameterHashByType diameters=new DiameterHashByType();
    diameters.setDiameter(typeCore, params.coreSigma);
    displayBox0.setDiameterHash(diameters);
            // Color scheme with different color for each arm and orange core
            ColorScheme colorScheme = new ColorScheme() {
                // Define colors for arms (red, blue, green, purple, yellow, cyan, magenta, etc.)
                private final Color[] armColors = new Color[]{
                        new Color(231, 76, 60),   // Red
                        new Color(52, 152, 219),  // Blue
                        new Color(46, 204, 113),  // Green
                        new Color(155, 89, 182),  // Purple
                        new Color(241, 196, 15),  // Yellow
                        new Color(26, 188, 156),  // Cyan
                        new Color(230, 126, 34),  // Orange
                        new Color(192, 57, 43),   // Dark Red
                        new Color(142, 68, 173),  // Dark Purple
                        new Color(39, 174, 96),   // Dark Green
                        new Color(41, 128, 185),  // Dark Blue
                        new Color(243, 156, 18),  // Dark Yellow
                        new Color(211, 84, 0),    // Dark Orange
                        new Color(231, 76, 60),   // Red (repeat)
                        new Color(52, 152, 219),  // Blue (repeat)
                        new Color(46, 204, 113),  // Green (repeat)
                        new Color(155, 89, 182),  // Purple (repeat)
                        new Color(241, 196, 15),  // Yellow (repeat)
                        new Color(26, 188, 156),  // Cyan (repeat)
                        new Color(230, 126, 34),  // Orange (repeat)
                };

                private final Color coreColor = new Color(255, 215, 0); // Gold core

                @Override
                public Color getAtomColor(IAtom a) {
                    // Check if atom is the core (atom 0)
                    if (a.getType() == typeCore) {
                        return coreColor;
                    }

                    // For arm beads: calculate which arm based on atom index
                    // Structure: atom 0 = core, atoms 1 to armLength = arm 0,
                    //            atoms armLength+1 to 2*armLength = arm 1, etc.
                    int atomIndex = a.getLeafIndex();

                    // Calculate arm number (subtract 1 for core, divide by armLength)
                    int armIndex = (atomIndex - 1) / armLength;

                    return armColors[armIndex % armColors.length];
                }
            };
            displayBox0.setColorScheme(colorScheme);



    simGraphic.makeAndDisplayFrame();

    final DisplayTextBox stepsBox = new DisplayTextBox();
    stepsBox.setLabel("Steps");
    DataSourceCountSteps stepCounter = new DataSourceCountSteps(integrator);
    DataPumpListener stepPump = new DataPumpListener(stepCounter, stepsBox, 1000);
    integrator.getEventManager().addListener(stepPump);
    simGraphic.add(stepsBox);

    return;
}

        long t1 = System.nanoTime();

        // if running interactively, don't use the file
        ActivityIntegrate ai = new ActivityIntegrate(integrator, steps/10);
        sim.getController().runActivityBlocking(ai);

        long t2 = System.nanoTime();
        System.out.println("equilibration finished: "+(t2-t1)/1e9);
        System.out.println("Angle move step size    " + angleMove.getStepSize());
        if (shuffleMove!=null) System.out.println("Shuffle move step size    "+shuffleMove.getStepSize());

        integrator.getMoveManager().setEquilibrating(false);
        MeterRadiusGyration meterRg = new MeterRadiusGyration(sim.box());
        AccumulatorAverageFixed accRg = new AccumulatorAverageFixed(steps/1000);
        DataPumpListener pumpRg = new DataPumpListener(meterRg, accRg, 10);
        integrator.getEventManager().addListener(pumpRg);
        XYZWriter xyzWriter = new XYZWriter(sim.box());
        xyzWriter.setFileName("star.xyz");
        xyzWriter.setIsAppend(true);
        integrator.getEventManager().addListener(new IntegratorListenerAction(xyzWriter,(int)steps/100));
        ai = new ActivityIntegrate(integrator, steps);
        sim.getController().runActivityBlocking(ai);

        System.out.println();

        System.out.println("Angle move acceptance " + angleMove.getTracker().acceptanceProbability());
        if (shuffleMove!=null) System.out.println("Shuffle move acceptance " + shuffleMove.getTracker().acceptanceProbability());
        System.out.println();

        double avgRg = accRg.getData(accRg.AVERAGE).getValue(0);
        double errRg = accRg.getData(accRg.ERROR).getValue(0);
        double corRg = accRg.getData(accRg.BLOCK_CORRELATION).getValue(0);

        System.out.println("Rg2: "+avgRg+"   err: "+errRg+"  cor: "+corRg);
        System.out.println("Core sigma: " + params.coreSigma + "  Core epsilon: " + params.epsCore);
        System.out.println("Rg: " + Math.sqrt(avgRg) + "   err: " + errRg/(2*Math.sqrt(avgRg)));
        //new add
        long t3 = System.nanoTime();
        // NOW write to CSV (file I/O shouldn't be included in simulation timing)
        String outputFile = System.getProperty("user.dir") + "/results.csv";
        boolean fileExists = new java.io.File(outputFile).exists();

        try (FileWriter fw = new FileWriter(outputFile, true)) {
            if (!fileExists) {
                fw.write("numArms,armLength,coreSigma,epsCore,ideal,temperature,Rg2,Rg2_err,Rg,Rg_err,correlation,numSteps\n");
            }

            String idealStr = params.ideal ? "IDEAL" : String.format("%.2f", params.temperature);
            fw.write(String.format("%d,%d,%.2f,%.2f,%s,%s,%.8f,%.8f,%.8f,%.8f,%.8f,%d%n",
                    params.numArms,
                    params.armLength,
                    params.coreSigma,
                    params.epsCore,
                    params.ideal,
                    idealStr,
                    avgRg,
                    errRg,
                    Math.sqrt(avgRg),
                    errRg/(2*Math.sqrt(avgRg)),
                    corRg,
                    params.numSteps));
            System.out.println("Results written to: " + outputFile);
        } catch (Exception e) {
            System.err.println("Warning: Could not write to results.csv: " + e.getMessage());
        }

        System.out.println("time: "+(t3-t2)/1e9);
    }

    /**
     * Inner class for parameters
     */
    public static class VirialStarParams extends ParameterBase {
        public int armLength = 3;
        public double temperature = 1;
        public long numSteps = 1000000;
        public int numArms = 2;
        public boolean ideal = false;
        // NEW
        public double coreSigma = 1.0;  // σ_core (also used in LJ)
        public double epsCore   = 1.0;  // ε_core

        // HS reference diameter options (if needed later)
        public boolean useCoreInHSRef = true;
        public double hsExtra = 0.0;
    }
}
