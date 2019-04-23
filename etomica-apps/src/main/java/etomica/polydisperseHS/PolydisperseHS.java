/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.polydisperseHS;


import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;
import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHash;
import etomica.atom.IAtom;
import etomica.atom.iterator.AtomIteratorLeafFilteredType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.history.HistoryScrolling;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPressureHard;
import etomica.data.meter.MeterRDF;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorScheme;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorListenerAction;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.random.IRandom;
import etomica.zeolite.MSDCoordWriter;

import java.awt.*;

/**
 * Three-dimensional hard-sphere molecular dynamics simulation, using
 * neighbor listing.
 * <p>
 * Developed as a prototype and example for the construction of a basic simulation.
 *
 * @author David Kofke and Andrew Schultz
 */
public class PolydisperseHS extends Simulation {

    //the following fields are made accessible for convenience to permit simple
    //mutation of the default behavior

    /**
     * The Box holding the atoms.
     */
    public final Box box;
    /**
     * The Integrator performing the dynamics.
     */
    public final IntegratorHard integrator;
    /**
     * The single hard-sphere species.
     */
    public final SpeciesSpheresMono species;
    /**
     * The hard-sphere potential governing the interactions.
     */
    public final P2HardSpherePoly potential;

    public final PotentialMaster potentialMaster;

    public final ActivityIntegrate activityIntegrate;

    /**
     * Makes a simulation according to the specified parameters.
     * @param params Parameters as defined by the inner class HSMD3DParam
     */
    public PolydisperseHS(HSPolyParam params) {

        super(Space3D.getInstance());

        species = new SpeciesSpheresMono(this, space);
        species.setIsDynamic(true);
        addSpecies(species);

        box = this.makeBox();

        double neighborRangeFac = 1.6;

        double [] sigmas = new double[params.numAtoms];
        double vsum = 0;

        if(params.isPolydisperse){
            for(int i=0;i<params.numAtoms;i++) {
                double r = random.nextGaussian();
                sigmas[i] = 1 + params.polySigma * r;
                vsum += Math.PI * sigmas[i]*sigmas[i]*sigmas[i]/6;
            }
        }else{//binary
            int [] randAtomList = new int[params.numAtoms];
            for(int i=0;i<params.numAtoms;i++){
                randAtomList[i] = i;
            }
            for(int i=0;i<params.numAtoms-1;i++){//0,1,...N-2
                int iRand =random.nextInt(params.numAtoms-i); // iRand = 0,1,...,(N-1)-i.
                int t = randAtomList[params.numAtoms-1-i];
                randAtomList[params.numAtoms-1-i] = randAtomList[iRand];
                randAtomList[iRand] = t;
            }
            int numAtomsA = (int)Math.round(params.binaryComp*params.numAtoms);

            for(int i=0;i<params.numAtoms;i++) {
                int iRand = randAtomList[i];
                if(i<numAtomsA){
                    sigmas[iRand] = params.sigmaRatio; // sigmaA
                }else{
                    sigmas[iRand] = 1.0; // sigmaB
                }
                vsum += Math.PI * sigmas[iRand]*sigmas[iRand]*sigmas[iRand]/6;
            }
        }
        double vBox = vsum / params.initEta;
        potential = new P2HardSpherePoly(space, sigmas,true);

        boolean useNeighborLists = true;

        potentialMaster = useNeighborLists ? new PotentialMasterList(this, potential.getRange() * neighborRangeFac, space) : new PotentialMasterMonatomic(this);

        integrator = new IntegratorHard(this, potentialMaster, box);
        integrator.setIsothermal(false);
        integrator.setTimeStep(params.timestep);



        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);

        AtomType leafType = species.getLeafType();


        potentialMaster.addPotential(potential, new AtomType[]{leafType, leafType});

        box.setNMolecules(species, params.numAtoms);
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(params.numAtoms / vBox);
        inflater.actionPerformed();

        System.out.println(box.getNMolecules(species)+" atoms");
        System.out.println(box.getBoundary().getEdgeVector(0).getX(0)+" "+box.getBoundary().getEdgeVector(1).getX(1)+" "+box.getBoundary().getEdgeVector(2).getX(2));

        if (space.D() == 3) {
            new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);
        } else {
            new ConfigurationLattice(new LatticeOrthorhombicHexagonal(space), space).initializeCoordinates(box);
        }

        if (useNeighborLists) {
            NeighborListManager nbrManager = ((PotentialMasterList) potentialMaster).getNeighborManager(box);
            integrator.getEventManager().addListener(nbrManager);
        } else {
            integrator.getEventManager().addListener(new IntegratorListenerAction(new BoxImposePbc(box, space)));
        }
        integrator.reset();
    }




    public static void main(String[] args) {
        final String APP_NAME = "PolyHS";

        HSPolyParam params = new HSPolyParam();
        final PolydisperseHS sim = new PolydisperseHS(params);

        MeterPressureHard meterP = new MeterPressureHard(sim.integrator);
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.potentialMaster, sim.box);

        // RDF g(r)
        MeterRDF meterRDF = new  MeterRDF(sim.space);
        meterRDF.setBox(sim.box);
        meterRDF.getXDataSource().setXMax(params.rMaxRDF);
        meterRDF.getXDataSource().setNValues(params.nBinsRDF);



        IAction init = new IAction() {
            public void actionPerformed() {
                if (params.initEta >= params.finalEta) return;

                double prevEta = params.initEta;
                double deltaEta = 0.01;
                for (double eta = params.initEta+deltaEta; eta<params.finalEta+deltaEta; eta+=deltaEta) {
                    if (eta > params.finalEta) eta = params.finalEta;
                    System.out.println(prevEta +"(current) => "+eta);
                    sim.activityIntegrate.setMaxSteps(params.numStepsComp);
                    sim.activityIntegrate.actionPerformed();

                    System.out.println("uncompressed density: " + params.numAtoms / sim.box.getBoundary().volume());
                    System.out.println(" uncompressed energy: " + (Math.log(meterPE.getDataAsScalar())));

                    BoxInflate inflater = new BoxInflate(sim.box, sim.space);
                    inflater.setScale(Math.pow(eta/ prevEta, -1.0 / 3.0));
                    inflater.actionPerformed();

                    sim.integrator.reset();

                    System.out.println("compressed density : " + params.numAtoms / sim.box.getBoundary().volume());
                    System.out.println("  compressed energy: " + (-Math.log(meterPE.getDataAsScalar())));

                    System.out.println("---------------------------");

                    sim.integrator.resetStepCount();
                    meterP.reset();

                    prevEta = eta;

                    if (eta==params.finalEta){ //for Graphics
                        System.out.println(prevEta+"(final)");
                        sim.activityIntegrate.setMaxSteps(Long.MAX_VALUE);
                        break; // then go to scond line below: sim.getController().addAction(sim.activityIntegrate);
                    }

                }
            }
        };




        if(params.isGraphics){
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME);
            simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));
            simGraphic.makeAndDisplayFrame(APP_NAME);
            DiameterHash diameterHash = new DiameterHash() {
                @Override
                public double getDiameter(IAtom atom) {
                    return sim.potential.getCollisionDiameter()[atom.getLeafIndex()];
                }
            };
            simGraphic.getDisplayBox(sim.box).setDiameterHash(diameterHash);
            ColorScheme colorScheme = new ColorScheme() {
                protected Color[] allColors;
                public Color getAtomColor(IAtom a) {
                    if (allColors==null) {
                        allColors = new Color[768];
                        for (int i=0; i<256; i++) {
                            allColors[i] = new Color(255-i,i,0);
                        }
                        for (int i=0; i<256; i++) {
                            allColors[i+256] = new Color(0,255-i,i);
                        }
                        for (int i=0; i<256; i++) {
                            allColors[i+512] = new Color(i,0,255-i);
                        }
                    }
                    return allColors[(2*a.getLeafIndex()) % 768];
                }
            };
            simGraphic.getDisplayBox(sim.box).setColorScheme(colorScheme);

            sim.getController().removeAction(sim.activityIntegrate);
            sim.getController().addAction(init);
            sim.getController().addAction(sim.activityIntegrate);

            AccumulatorHistory historyP = new AccumulatorHistory(new HistoryCollapsingAverage());
            DataPumpListener pumpP = new DataPumpListener(meterP, historyP, params.numAtoms);
            sim.integrator.getEventManager().addListener(pumpP);
            DisplayPlot plotP = new DisplayPlot();
            historyP.setDataSink(plotP.getDataSet().makeDataSink());
            plotP.setLabel("P");
            simGraphic.add(plotP);


            AccumulatorHistory historyPE = new AccumulatorHistory(new HistoryScrolling());
            DataPumpListener pumpPE = new DataPumpListener(meterPE, historyPE, params.numAtoms);
            sim.integrator.getEventManager().addListener(pumpPE);
            DisplayPlot plotPE = new DisplayPlot();
            historyPE.setDataSink(plotPE.getDataSet().makeDataSink());
            plotPE.setLabel("PE");
            simGraphic.add(plotPE);

            // RDF
            DisplayPlot plotRDF = new DisplayPlot();
            DataPumpListener pumpRDF = new DataPumpListener(meterRDF,plotRDF.getDataSet().makeDataSink(),10);
            sim.integrator.getEventManager().addListener(pumpRDF);
            sim.integrator.getEventManager().addListener(new IntegratorListenerAction(meterRDF));
            plotRDF.setLabel("RDF");
            simGraphic.add(plotRDF);

        }else{
        //Compress from initEta to finalEta
            init.actionPerformed();
        //Equilibarate with finalEta
            sim.activityIntegrate.setMaxSteps(params.numStepsEqu);
            sim.activityIntegrate.actionPerformed();

            System.out.println("=========================================================================");
            System.out.println("After Equilibaration: ");
            System.out.println("density: " + params.numAtoms / sim.box.getBoundary().volume());
            System.out.println(" energy: " + (Math.log(meterPE.getDataAsScalar())));
            System.out.println("=========================================================================");

            sim.integrator.resetStepCount();
            meterP.reset();


        //Production

        // coorWriter
            String filename_pos = "coordinates_uw";
            System.out.println(filename_pos);
            MSDCoordWriter coordWriter = new MSDCoordWriter(sim.integrator, sim.box, filename_pos, params.writeIntervalPos);
            coordWriter.setIterator(new AtomIteratorLeafFilteredType(sim.box, sim.species.getLeafType()));
            System.out.println("created MSDCoordWriter");
            sim.getController().getEventManager().addListener(coordWriter);

        //velWriter
            String filename_vel = "velocities";
            System.out.println(filename_vel);
            VelocityWriter velWriter = new VelocityWriter(sim.integrator, sim.box, filename_vel, params.writeIntervalVel);
            velWriter.setIterator(new AtomIteratorLeafFilteredType(sim.box, sim.species.getLeafType()));
            System.out.println("created velCoordWriter");
            sim.getController().getEventManager().addListener(velWriter);


            int numBlocks = 100;
            int interval = 10;
            long blockSize = params.numStepsProd/(numBlocks*interval);
            if (blockSize == 0) blockSize = 1;
            System.out.println("block size "+blockSize+" interval_stat "+interval);
            AccumulatorAverageFixed accumulatorP = new AccumulatorAverageFixed(blockSize);
            DataPumpListener accumulatorPPump = new DataPumpListener(meterP, accumulatorP, interval);
            sim.integrator.getEventManager().addListener(accumulatorPPump);

            // RDF
            sim.integrator.getEventManager().addListener(new IntegratorListenerAction(meterRDF));


        //Run ...
            sim.activityIntegrate.setMaxSteps(params.numStepsProd);
            sim.getController().actionPerformed();

            // statistics
            DataGroup dataP = (DataGroup)accumulatorP.getData();
            IData dataPAvg = dataP.getData(accumulatorP.AVERAGE.index);
            IData dataPErr = dataP.getData(accumulatorP.ERROR.index);
            IData dataPCorrelation = dataP.getData(accumulatorP.BLOCK_CORRELATION.index);
            double pAvg = dataPAvg.getValue(0);
            double pErr = dataPErr.getValue(0);
            double pCor = dataPCorrelation.getValue(0);
            System.out.println(" P = "+pAvg+"  Err = "+pErr+"  cor: "+pCor);


//            IData dataRDF = meterRDF.getData();
//            IData dataR = meterRDF.getIndependentData(0);
//            for (int i = 0; i < dataRDF.getLength(); i++){
//                double r = dataR.getValue(i);
//                double g = dataRDF.getValue(i);
//                System.out.println(r + " "+ g);
//            }


        }

    }

    public static class HSPolyParam extends ParameterBase {
        public boolean isGraphics = true;
        public int nC = 4;
        public int numAtoms = 4*nC*nC*nC;

        //Compression
        public double initEta    = 0.55;
        public double finalEta   = 0.57;

        //Polydisperse
        public boolean isPolydisperse = !true;


        public double polySigma = 0.1;


        public double sigmaRatio = 1.1;
        public double binaryComp = 0.2;

        public long numStepsComp = 10000;
        public long numStepsEqu  = 10000;

        public double timestep   = 0.01;
        public long numStepsProd = 100000;

        //RDF
        public double rMaxRDF = 5.0;
        public int nBinsRDF = 100;

        //Coordinates
        public int writeIntervalPos = 10;
        //Velocites
        public int writeIntervalVel = 10;

    }
}