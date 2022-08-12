/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.config.ConformationGeneric;
import etomica.config.IConformation;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountSteps;
import etomica.data.IData;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorScheme;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.potential.P1Anharmonic;
import etomica.potential.P2Harmonic;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.compute.PotentialCompute;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputeField;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space1d.Space1D;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
public class SimQuantumAO extends Simulation {

    public SpeciesGeneral species;
    public Box box;
    public IntegratorMC integrator;
    public PotentialComputeField pcP1;
    public PotentialMasterBonding pmBonding;
    public PotentialCompute pm;
    public double betaN;
    public double mass;
    public double k2_kin;

    public static final double kB = 1.0;//Constants.BOLTZMANN_K



    public static final double hbar = 1.0; //Constants.PLANCK_H/(2.0*Math.PI);





    public SimQuantumAO(Space space, int nBeads, double temperature, double omega, double k4) {
        super(space);
        Vector[] initCoords = new Vector[nBeads];
        for (int i = 0; i < nBeads; i++){
            initCoords[i] = space.makeVector();
        }
        IConformation conformation = new ConformationGeneric(initCoords);

        SpeciesGeneral species = new SpeciesBuilder(space)
                .addCount(AtomType.simpleFromSim(this), nBeads)
                .withConformation(conformation)
                .build();
        addSpecies(species);

        box = this.makeBox(new BoundaryRectangularNonperiodic(space));
        box.setNMolecules(species, 1);
        box.getBoundary().setBoxSize(Vector.of(new double[]{1}));
        mass = species.getLeafType().getMass();

        //pm2 that uses the full PI potential, for data collection
        //spring P2 part (x_i-x_{i+1})^2
        pmBonding = new PotentialMasterBonding(getSpeciesManager(), box);

        double beta = 1.0/(kB*temperature);
        betaN = beta/nBeads;
        double omegaN = 1.0/(hbar*betaN);

        k2_kin = mass*omegaN*omegaN/nBeads;

        P2Harmonic p2Bond = new P2Harmonic(k2_kin, 0);
        List<int[]> pairs = new ArrayList<>();
        for (int i=0; i<nBeads; i++) {
            int[] p = new int[]{i,i+1};
            for (int j=0; j<p.length; j++) {
                if (p[j] >= nBeads) p[j] -= nBeads;
            }
            pairs.add(p);
        }
        pmBonding.setBondingPotentialPair(species, p2Bond, pairs);

        //potential p1 part
        double k2 = mass*omega*omega;
        P1Anharmonic p1ah = new P1Anharmonic(space, k2/nBeads, k4/nBeads);
        pcP1 = new PotentialComputeField(getSpeciesManager(), box);
        pcP1.setFieldPotential(species.getLeafType(), p1ah);

        //TOTAL
        pm = new PotentialComputeAggregate(pmBonding, pcP1);


        integrator = new IntegratorMC(pm, random, temperature, box);
        MCMoveHO  atomMove = new MCMoveHO(space, pm, random, temperature, omega, box);
//        MCMoveAtom atomMove = new MCMoveAtom(random, pm, box);
        integrator.getMoveManager().addMCMove(atomMove);




        getController().addActivity(new ActivityIntegrate(integrator));



        // make an IntegratorMC and write a move with a doTrial that samples the harmonic part of your system


    }

    public static void main(String[] args) {

        OctaneParams params = new OctaneParams();
        ParseArgs.doParseArgs(params, args);

        double temperature = params.temperature;
        double omega = params.omega;;
        double k4 = params.k4;;
        int nBeads = params.nBeads;
        boolean graphics = params.graphics;
        long numSteps = params.numSteps;

        final SimQuantumAO sim = new SimQuantumAO(Space1D.getInstance(), nBeads, temperature, omega, k4);
        sim.integrator.reset();

        System.out.println(numSteps+" steps");
        System.out.println(" nBeads: " + nBeads);
        System.out.println(" temperature: " + temperature);
        System.out.println(" mass: " + sim.mass + " omega: " + omega + " k4: " + k4);
        System.out.println(" k2_kin: " + sim.k2_kin);

        MeterMSDHO meterMSDHO = new MeterMSDHO(nBeads, sim.box);
        MeterPIPrim meterPrim = new MeterPIPrim(sim.pmBonding, sim.pcP1, sim.betaN, nBeads);
        MeterPIVir meterVir = new MeterPIVir(sim.pcP1, sim.betaN, nBeads, sim.box);
        MeterPIHMA meterHMA = new MeterPIHMA(sim.pmBonding, sim.pcP1, sim.betaN, nBeads, omega, sim.box);





        if (graphics) {
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            simGraphic.setPaintInterval(sim.box, 1000);
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

            DisplayTextBox timer = new DisplayTextBox();
            DataSourceCountSteps counter = new DataSourceCountSteps(sim.integrator);
            DataPumpListener counterPump = new DataPumpListener(counter, timer, 100);
            sim.integrator.getEventManager().addListener(counterPump);
            simGraphic.getPanel().controlPanel.add(timer.graphic());

            simGraphic.makeAndDisplayFrame(" PIMC ");

            return;
        }











        System.out.flush();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps/5));
        sim.integrator.getMoveManager().setEquilibrating(false);



        int numBlocks = 100;


        ///////////////////////////////////
        int interval = 10;
        ////////////////////////////////////////////////



        long blockSize = numSteps/numBlocks;
        if (blockSize == 0) blockSize = 1;
        System.out.println("block size "+blockSize);

        AccumulatorAverageFixed accumulatorMSD = new AccumulatorAverageFixed(blockSize);
        DataPumpListener accumulatorPumpMSD = new DataPumpListener(meterMSDHO, accumulatorMSD, interval);
        sim.integrator.getEventManager().addListener(accumulatorPumpMSD);

        // Primitive E_n
        AccumulatorAverageFixed accumulatorPrim = new AccumulatorAverageFixed(blockSize);
        DataPumpListener accumulatorPumpPrim = new DataPumpListener(meterPrim, accumulatorPrim, interval);
        sim.integrator.getEventManager().addListener(accumulatorPumpPrim);

        // Virial E_n
        AccumulatorAverageFixed accumulatorVir = new AccumulatorAverageFixed(blockSize);
        DataPumpListener accumulatorPumpVir = new DataPumpListener(meterVir, accumulatorVir, interval);
        sim.integrator.getEventManager().addListener(accumulatorPumpVir);

        // HMA E_n
        AccumulatorAverageFixed accumulatorHMA = new AccumulatorAverageFixed(blockSize);
        DataPumpListener accumulatorPumpHMA = new DataPumpListener(meterHMA, accumulatorHMA, interval);
        sim.integrator.getEventManager().addListener(accumulatorPumpHMA);



        final long startTime = System.currentTimeMillis();

        //run
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps));


        DataGroup dataMSD = (DataGroup)accumulatorMSD.getData();
        IData dataMSDErr = dataMSD.getData(accumulatorMSD.ERROR.index);
        IData dataMSDAvg = dataMSD.getData(accumulatorMSD.AVERAGE.index);
        IData dataMSDCorrelation = dataMSD.getData(accumulatorMSD.BLOCK_CORRELATION.index);
        double avgMSD = dataMSDAvg.getValue(0);
        double errMSD = dataMSDErr.getValue(0);
        double corMSD = dataMSDCorrelation.getValue(0);



        DataGroup dataPrim = (DataGroup)accumulatorPrim.getData();
        IData dataErrPrim = dataPrim.getData(accumulatorPrim.ERROR.index);
        IData dataAvgPrim = dataPrim.getData(accumulatorPrim.AVERAGE.index);
        IData dataCorrelationPrim = dataPrim.getData(accumulatorPrim.BLOCK_CORRELATION.index);
        double avgEnPrim = dataAvgPrim.getValue(0);
        double errEnPrim = dataErrPrim.getValue(0);
        double corEnPrim = dataCorrelationPrim.getValue(0);
        System.out.println("\nEn_primitive: " + avgEnPrim  + " +/- " + errEnPrim + " cor: " + corEnPrim);


        DataGroup dataVir = (DataGroup)accumulatorVir.getData();
        IData dataErrVir = dataVir.getData(accumulatorVir.ERROR.index);
        IData dataAvgVir = dataVir.getData(accumulatorVir.AVERAGE.index);
        IData dataCorrelationVir = dataVir.getData(accumulatorVir.BLOCK_CORRELATION.index);
        double avgEnVir = dataAvgVir.getValue(0);
        double errEnVir = dataErrVir.getValue(0);
        double corEnVir = dataCorrelationVir.getValue(0);
        System.out.println("En_virial: " + avgEnVir  + " +/- " + errEnVir + " cor: " + corEnVir);


        DataGroup dataHMA = (DataGroup)accumulatorHMA.getData();
        IData dataErrHMA = dataHMA.getData(accumulatorHMA.ERROR.index);
        IData dataAvgHMA = dataHMA.getData(accumulatorHMA.AVERAGE.index);
        IData dataCorrelationHMA = dataHMA.getData(accumulatorHMA.BLOCK_CORRELATION.index);
        double avgEnHMA = dataAvgHMA.getValue(0);
        double errEnHMA = dataErrHMA.getValue(0);
        double corEnHMA = dataCorrelationHMA.getValue(0);
        System.out.println("En_HMA: " + avgEnHMA  + " +/- " + errEnHMA + " cor: " + corEnHMA);




        System.out.println("\n Theory");
        System.out.println("MSD_sim: " + avgMSD + " +/- " + errMSD + " cor: " + corMSD);
        System.out.println("MSDc: " + sim.kB*temperature/ sim.mass/omega/omega);

        System.out.println("MSDq: " + hbar/sim.mass/omega*(0.5+1.0/(Math.exp(hbar*omega/temperature)-1.0)));

        double EnQ = hbar*omega*(0.5 + 1/(Math.exp(nBeads*sim.betaN*hbar*omega)-1.0));
        double EnC = temperature;
        System.out.println("\nEnQ: " + EnQ + " EnC: " + EnC);



        long endTime = System.currentTimeMillis();
        System.out.println();
        System.out.println("time: " + (endTime - startTime)/1000.0);

    }

    public static class OctaneParams extends ParameterBase {
        public double temperature = 1.0;
        public int nBeads = 111; //must be odd for now!
        public boolean graphics = false;
        public double omega = 1.0; // k2=m*w^2
        public double k4 = 24.0;
        public long numSteps = 1_000_000;
    }
}