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
import etomica.data.IData;


import etomica.data.meter.MeterRadiusGyrationComponents;
import etomica.data.types.DataVector;
import etomica.graphics.*;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.potential.*;
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
import etomica.virial.mcmove.MCMoveClusterAngleMulti;
import etomica.virial.mcmove.MCMoveClusterStretch;
import etomica.virial.simulations.theta.MCMoveClusterShuffle;

import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import etomica.zeno.MeterZENO;
import etomica.zeno.WalkerExterior;



public class VirialStarSingle {

    public static void main(String[] args) {
        VirialStarParams params = new VirialStarParams();

       if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.armLength = 60;
            params.numArms = 4;
            params.temperature = 10;
            params.numSteps = 200000;
            params.coreSigma = 1.0;  //new
            params.epsCore = 1.0;// new
            params.ideal = false;
            params.rc = 2.5;
            params.zenoWalks=5000;

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
        TruncationFactory tf = new TruncationFactorySimple(params.rc);

        //IPotential2 p2CoreCore = new P2LennardJones(epsC,  sigC);//new
        //IPotential2 p2MonoMono = new P2LennardJones(epsM,  sigM);//new
       // IPotential2 p2CoreMono = new P2LennardJones(epsCM, sigCM);//new

        double coreRadius = params.coreSigma / 2.0;

        IPotential2 p2CoreCore = new P2ParticleParticle(coreRadius, params.sigmaAtom, params.epsilonAtom, params.density);
        IPotential2 p2MonoMono = P2LennardJones.makeTruncated(sigM, epsM,tf);
        IPotential2 p2CoreMono = new P2ParticleSegment(coreRadius, params.sigmaAtom, params.epsilonAtom, params.density);

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
        MCMoveClusterAngle angleMoveInner=null;
        if (armLength>2){
             angleMoveInner = new MCMoveClusterAngle(pc, space, bonding, sim.getRandom(), 1,
                new MCMoveClusterAngle.AtomChooserStarfl(sim.getRandom(), numArms, armLength,2, (armLength+1)/2));
        angleMoveInner.setBox(sim.box());
        integrator.getMoveManager().addMCMove(angleMoveInner,0.1);}
        MCMoveClusterAngle angleMoveOuter=new MCMoveClusterAngle(pc, space, bonding, sim.getRandom(), 1,
                new MCMoveClusterAngle.AtomChooserStarfl(sim.getRandom(), numArms, armLength, (armLength+1) / 2 + 1, armLength));

        angleMoveOuter.setBox(sim.box());
        angleMoveOuter.setFixedCOM(false);
        integrator.getMoveManager().addMCMove(angleMoveOuter, 0.1);
        MCMoveClusterAngleMulti angleMoveMultiInner=null;
        if (armLength>2) {
            angleMoveMultiInner = new MCMoveClusterAngleMulti(pc, space, bonding, sim.getRandom(), 1,
                    new MCMoveClusterAngle.AtomChooserStarfl(sim.getRandom(), numArms, armLength, 2, (armLength+1) / 2), 5);

            angleMoveMultiInner.setBox(sim.box());

            angleMoveMultiInner.setFixedCOM(false);

            integrator.getMoveManager().addMCMove(angleMoveMultiInner, 0.4);
        }

        MCMoveClusterAngleMulti angleMoveMultiOuter = new MCMoveClusterAngleMulti(pc, space, bonding, sim.getRandom(), 1,
                new MCMoveClusterAngle.AtomChooserStarfl(sim.getRandom(), numArms, armLength, (armLength+1) / 2 + 1, armLength), 10);
        angleMoveMultiOuter.setBox(sim.box());

        angleMoveMultiOuter.setFixedCOM(false);

        integrator.getMoveManager().addMCMove(angleMoveMultiOuter, 0.4);



        /*MCMoveClusterShuffle shuffleMove = null;
        if (armLength > 6 ) {
            shuffleMove = new MCMoveClusterShuffle(pc, space, sim.getRandom());
            shuffleMove.setExcludedAtom(0);
            shuffleMove.setBox(sim.box());
            shuffleMove.setBonding(bonding);
            shuffleMove.setStepSizeMax(armLength - 2);
            integrator.getMoveManager().addMCMove(shuffleMove);
            ((MCMoveStepTracker) shuffleMove.getTracker()).setAcceptanceTarget(0.3);

        }*/


        if (false) {
    //double size = (2*armLength + 5) * 1.5;
    double size = 50;
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
        if(angleMoveInner!=null) System.out.println("Angle move step size Inner   " + angleMoveInner.getStepSize());
        System.out.println("Angle move step size  Outer  " + angleMoveOuter.getStepSize());
        if(angleMoveMultiInner!=null) System.out.println("Angle move step size MultiInner   " + angleMoveMultiInner.getStepSize());
        System.out.println("Angle move step size MultiOuter   " + angleMoveMultiOuter.getStepSize());
       // if (shuffleMove!=null) System.out.println("Shuffle move step size    "+shuffleMove.getStepSize());

        integrator.getMoveManager().setEquilibrating(false);
        MeterRadiusGyration meterRg = new MeterRadiusGyration(sim.box());
        MeterRadiusGyrationComponents componentsRg= new MeterRadiusGyrationComponents(sim.box());
        AccumulatorAverageFixed accRg = new AccumulatorAverageFixed(steps/1000);
        AccumulatorAverageFixed accComponents = new AccumulatorAverageFixed(steps/1000);
        DataPumpListener pumpRg = new DataPumpListener(meterRg, accRg, 10);
        DataPumpListener pumpcomponentsRg = new DataPumpListener(componentsRg,accComponents,10);
        //new
        double coreWeight = params.density * (4.0/3.0) * Math.PI * Math.pow(coreRadius, 3);
        MeterRadiusGyrationComponentsCore componentsRgCore =
                new MeterRadiusGyrationComponentsCore(sim.box(), typeCore, coreRadius, coreWeight);
        AccumulatorAverageFixed accComponentsCore = new AccumulatorAverageFixed(steps/1000);
        DataPumpListener pumpComponentsRgCore = new DataPumpListener(componentsRgCore, accComponentsCore, 10);

        MeterArmProfile meterProfile = new MeterArmProfile(sim.box(), numArms, armLength);
        AccumulatorAverageFixed accProfile = new AccumulatorAverageFixed(steps / 1000);
        DataPumpListener pumpProfile = new DataPumpListener(meterProfile, accProfile, 10);
        integrator.getEventManager().addListener(pumpProfile);
        integrator.getEventManager().addListener(pumpRg);
        integrator.getEventManager().addListener(pumpcomponentsRg);
        integrator.getEventManager().addListener(pumpComponentsRgCore); //new for core


        // Hydrodynamic radius (Rh) using ZENO

        double[] zenoSigma = {params.coreSigma, 1.0};
        MeterZENO meterZeno =
                new MeterZENO(sim.box(), sim.getRandom(), zenoSigma, false, WalkerExterior.Shape.SPHERE);
        meterZeno.setEigenvaluesForPade(true);

        long zenoWalks = params.zenoWalks;
        int zenoInterval = params.zenoInterval;
         //long zenoWalks = 5000;
         //int zenoInterval = 1000;


        meterZeno.setNumWalks(zenoWalks);
// charge distribution: COM weights same as core-aware Rg
        double[] massByType = new double[2];
        massByType[typeCore.getIndex()] = coreWeight;
        massByType[typeMono.getIndex()] = 1.0;
        meterZeno.setMassByType(massByType);
        meterZeno.setChargeBinWidth(params.chargeBinWidth);
/*     MeterZENO.MeterHydrodynamicRadius meterRh =
                meterZeno.makeMeterHydrodynamicRadius();

        AccumulatorAverageFixed accRh =
                new AccumulatorAverageFixed(10);

        DataPumpListener pumpRh =
                new DataPumpListener(meterRh, accRh, zenoInterval);

        integrator.getEventManager().addListener(pumpRh); */
        AccumulatorAverageFixed accZeno =
                new AccumulatorAverageFixed(10);

        DataPumpListener pumpZeno =
                new DataPumpListener(meterZeno, accZeno, zenoInterval);

        integrator.getEventManager().addListener(pumpZeno);

        XYZWriter xyzWriter = new XYZWriter(sim.box());
        xyzWriter.setFileName("star.xyz");
        xyzWriter.setIsAppend(true);
        integrator.getEventManager().addListener(new IntegratorListenerAction(xyzWriter,(int)steps/100));
        ai = new ActivityIntegrate(integrator, steps);
        sim.getController().runActivityBlocking(ai);

        System.out.println();

        if(angleMoveInner!=null) System.out.println("Angle move acceptance Inner   " + angleMoveInner.getTracker().acceptanceProbability());
        System.out.println("Angle move acceptance  Outer  " + angleMoveOuter.getTracker().acceptanceProbability());
        if(angleMoveMultiInner!=null) System.out.println("Angle move acceptance MultiInner   " + angleMoveMultiInner.getTracker().acceptanceProbability());
        System.out.println("Angle move acceptance MultiOuter   " + angleMoveMultiOuter.getTracker().acceptanceProbability());


        //if (shuffleMove!=null) System.out.println("Shuffle move acceptance " + shuffleMove.getTracker().acceptanceProbability());
        System.out.println();

         double avgRg = accRg.getData(accRg.AVERAGE).getValue(0);
        double errRg = accRg.getData(accRg.ERROR).getValue(0);
        double corRg = accRg.getData(accRg.BLOCK_CORRELATION).getValue(0);
        IData zenoAvg = accZeno.getData(accZeno.AVERAGE);
        IData zenoErr = accZeno.getData(accZeno.ERROR);
        IData zenoCor = accZeno.getData(accZeno.BLOCK_CORRELATION);
        // index 0 = Rh, index 5 = intrinsic viscosity (per bead)
        double avgRh = zenoAvg.getValue(0);
        double errRh = zenoErr.getValue(0);
        double corRh = zenoCor.getValue(0);
        double avgEta = zenoAvg.getValue(5);
        double errEta = zenoErr.getValue(5);
        double corEta = zenoCor.getValue(5);

        DataVector componentsavgRg = (DataVector) accComponents.getData(accComponents.AVERAGE);
        DataVector componentserrRg = (DataVector) accComponents.getData(accComponents.ERROR);
        DataVector componentscorRg = (DataVector) accComponents.getData(accComponents.BLOCK_CORRELATION);
        //new for core

        IData avgCore = accComponentsCore.getData(accComponentsCore.AVERAGE);
        IData errCore = accComponentsCore.getData(accComponentsCore.ERROR);
        IData corCore = accComponentsCore.getData(accComponentsCore.BLOCK_CORRELATION);

        double Rg2xCore = avgCore.getValue(0);
        double Rg2yCore = avgCore.getValue(1);
        double Rg2zCore = avgCore.getValue(2);
        double Rg2Core  = avgCore.getValue(3);
        double Rg2xCore_err = errCore.getValue(0);
        double Rg2yCore_err = errCore.getValue(1);
        double Rg2zCore_err = errCore.getValue(2);
        double Rg2Core_err  = errCore.getValue(3);
        double Rg2xCore_cor = corCore.getValue(0);
        double Rg2yCore_cor = corCore.getValue(1);
        double Rg2zCore_cor = corCore.getValue(2);
        double Rg2Core_cor  = corCore.getValue(3);
        double RgCore       = Math.sqrt(Rg2Core);
        double RgCore_err   = Rg2Core_err / (2.0 * RgCore);


        //System.out.println("componentsRg2: "+componentsavgRg+" componentserr: "+componentserrRg+" componentscor: "+componentscorRg);
        double componentsRg2x= componentsavgRg.x.getX(0);
        double componentsRg2y= componentsavgRg.x.getX(1);
        double componentsRg2z= componentsavgRg.x.getX(2);
        double componentserrx= componentserrRg.x.getX(0);
        double componentserry= componentserrRg.x.getX(1);
        double componentserrz= componentserrRg.x.getX(2);
        double componentscorx= componentscorRg.x.getX(0);
        double componentscory= componentscorRg.x.getX(1);
        double componentscorz= componentscorRg.x.getX(2);

        System.out.println("componentsRg2x: "+componentsRg2x+" componentsRg2y: "+componentsRg2y+" componentsRg2z: "+componentsRg2z
                +" componentserrx: "+componentserrx+" componentserry: "+componentserry+" componentserrz: "+componentserrz
                +" componentscorx: "+componentscorx+" componentscory: "+componentscory+" componentscorz: "+componentscorz);



        System.out.println("Rg2: "+avgRg+"   err: "+errRg+"  cor: "+corRg);
        System.out.println("Rg2 (core-aware): " + Rg2Core
                + "   x: " + Rg2xCore + " y: " + Rg2yCore + " z: " + Rg2zCore);//new
        System.out.println("Core sigma: " + params.coreSigma + "  Core epsilon: " + params.epsCore);
        System.out.println("Rg: " + Math.sqrt(avgRg) + "   err: " + errRg/(2*Math.sqrt(avgRg)));
        System.out.println("Rg2 (core-aware): " + Rg2Core
                + "   x: " + Rg2xCore + " y: " + Rg2yCore + " z: " + Rg2zCore);//newcore

        // Hydrodynamic radius results

        /*double avgRh = accRh.getData(accRh.AVERAGE).getValue(0);
        double errRh = accRh.getData(accRh.ERROR).getValue(0);
        double corRh = accRh.getData(accRh.BLOCK_CORRELATION).getValue(0);*/

        System.out.println("Rh: " + avgRh
                + "   err: " + errRh
                + "   cor: " + corRh);
        System.out.println("Core type index: " + typeCore.getIndex());
        System.out.println("Mono type index: " + typeMono.getIndex());
        System.out.println("eta: " + avgEta + "   err: " + errEta + "   cor: " + corEta);


        IData avgProfile = accProfile.getData(accProfile.AVERAGE);
        IData errProfile = accProfile.getData(accProfile.ERROR);
        IData corProfile = accProfile.getData(accProfile.BLOCK_CORRELATION);

        System.out.println("Arm radial profile r(i): monomer index i (1=closest to core, armLength=tip)");
        System.out.println("i, r_avg, r_err, r_cor");
        for (int k = 1; k <= armLength; k++) {
            System.out.printf("%d, %.6f, %.6f, %.6f%n",
                    k,
                    avgProfile.getValue(k - 1),
                    errProfile.getValue(k - 1),
                    corProfile.getValue(k - 1));
        }
        //new add
        long t3 = System.nanoTime();
        // NOW write to CSV (file I/O shouldn't be included in simulation timing)
        String outputFile = System.getProperty("user.dir") + "/results.csv";
        boolean fileExists = new java.io.File(outputFile).exists();

        try (FileWriter fw = new FileWriter(outputFile, true)) {
            if (!fileExists) {
                fw.write("numArms,armLength,coreSigma,epsCore,ideal,temperature,Rg2,Rg2_err,Rg,Rg_err,correlation," +
                        "Rg2x,Rg2x_err,Rg2x_cor,Rg2y,Rg2y_err,Rg2y_cor,Rg2z,Rg2z_err,Rg2z_cor," +
                        "Rg2xCore,Rg2xCore_err,Rg2xCore_cor,Rg2yCore,Rg2yCore_err,Rg2yCore_cor," +
                        "Rg2zCore,Rg2zCore_err,Rg2zCore_cor,Rg2Core,Rg2Core_err,Rg2Core_cor,RgCore,RgCore_err,numSteps\n");
            }

            String idealStr = params.ideal ? "IDEAL" : String.format("%.2f", params.temperature);
            fw.write(String.format("%d,%d,%.2f,%.2f,%s,%s,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%d%n",
                    params.numArms, params.armLength, params.coreSigma, params.epsCore, params.ideal, idealStr,
                    avgRg, errRg, Math.sqrt(avgRg), errRg/(2*Math.sqrt(avgRg)), corRg,
                    componentsRg2x, componentserrx, componentscorx,
                    componentsRg2y, componentserry, componentscory,
                    componentsRg2z, componentserrz, componentscorz,
                    Rg2xCore, Rg2xCore_err, Rg2xCore_cor,
                    Rg2yCore, Rg2yCore_err, Rg2yCore_cor,
                    Rg2zCore, Rg2zCore_err, Rg2zCore_cor,
                    Rg2Core, Rg2Core_err, Rg2Core_cor, RgCore, RgCore_err,
                    params.numSteps));
            System.out.println("Results written to: " + outputFile);
        } catch (Exception e) {
            System.err.println("Warning: Could not write to results.csv: " + e.getMessage());
        }
        // ADDED: separate profile.csv — one row per arm bead index per run
        // columns: numArms, armLength, coreSigma, epsCore, temperature, i, r_avg, r_err, r_cor
        String profileFile = System.getProperty("user.dir") + "/profile.csv";
        boolean profileFileExists = new java.io.File(profileFile).exists();
        try (FileWriter fw = new FileWriter(profileFile, true)) {
            if (!profileFileExists) {
                fw.write("numArms,armLength,coreSigma,epsCore,temperature,i,r_avg,r_err,r_cor\n");
            }
            for (int k = 1; k <= armLength; k++) {
                fw.write(String.format("%d,%d,%.2f,%.2f,%.4f,%d,%.8f,%.8f,%.8f%n",
                        numArms,
                        armLength,
                        params.coreSigma,
                        params.epsCore,
                        temperature,
                        k,
                        avgProfile.getValue(k - 1),
                        errProfile.getValue(k - 1),
                        corProfile.getValue(k - 1)));
            }
            System.out.println("Profile written to: " + profileFile);
        } catch (Exception e) {
            System.err.println("Warning: Could not write to profile.csv: " + e.getMessage());
        }
        // --- charge distribution P(r) ---
        double[] bins = meterZeno.getChargeBins();
        long nConf = meterZeno.getChargeConfs();
        double dr = params.chargeBinWidth;
        double totalCharge = 0;
        String chargeFile = System.getProperty("user.dir") + "/charge.csv";
        boolean chargeExists = new java.io.File(chargeFile).exists();
        try (FileWriter fw = new FileWriter(chargeFile, true)) {
            if (!chargeExists) {
                fw.write("numArms,armLength,coreSigma,epsCore,ideal,temperature,Rh,binWidth,nConf,r,r_over_Rh,C3P\n");
            }
            for (int k = 0; k < bins.length; k++) {
                double r = (k + 0.5) * dr;                     // bin center
                double q = bins[k] / nConf;                    // avg charge in shell
                totalCharge += q;
                double rho = q / (4 * Math.PI * r * r * dr);   // charge density
                double y = avgRh * avgRh * rho;                // = C^3 P(r), P with unit total charge
                fw.write(String.format("%d,%d,%.2f,%.2f,%s,%.4f,%.8f,%.4f,%d,%.6f,%.8f,%.8f%n",
                        numArms, armLength, params.coreSigma, params.epsCore, params.ideal, temperature,
                        avgRh, dr, nConf, r, r / avgRh, y));
            }
            System.out.println("Charge distribution written to: " + chargeFile);
        } catch (Exception e) {
            System.err.println("Warning: Could not write to charge.csv: " + e.getMessage());
        }
        System.out.println("charge check: total charge " + totalCharge + " vs Rh " + avgRh);
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
        public long zenoWalks = 2000;
        public int zenoInterval = 1000;
        public double chargeBinWidth = 0.1;   // bin width (sigma units) for charge distribution
        // NEW
        public double coreSigma = 1.0;  // σ_core (also used in LJ)
        public double epsCore   = 1.0;  // ε_core
        public double rc=4;

        // HS reference diameter options (if needed later)
        public boolean useCoreInHSRef = true;
        public double hsExtra = 0.0;
        //new
        public double sigmaAtom = 1.0;
        public double epsilonAtom = 1.0;
        public double density = 1;
    }
}
