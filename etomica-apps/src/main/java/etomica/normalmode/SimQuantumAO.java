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
import etomica.util.Constants;
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
        double mass = nBeads*species.getMass();

        //pm2 that uses the full PI potential, for data collection
        //spring P2 part (x_i-x_{i+1})^2
        pmBonding = new PotentialMasterBonding(getSpeciesManager(), box);
        double hbar = Constants.PLANCK_H/(2*Math.PI);
        double beta = 1.0/(Constants.BOLTZMANN_K*temperature);
        betaN = beta/nBeads;
        double omegaN = 1.0/(hbar*betaN);

        double k2_kin = mass*omegaN*omegaN/nBeads;

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
        P1Anharmonic p1ah = new P1Anharmonic(space, k2, k4);
        pcP1 = new PotentialComputeField(getSpeciesManager(), box);
        pcP1.setFieldPotential(species.getLeafType(), p1ah);

        //TOTAL
        pm = new PotentialComputeAggregate(pmBonding, pcP1);


        integrator = new IntegratorMC(pm, random, temperature, box);
        MCMoveHO  atomMove = new MCMoveHO(space, pm, random, temperature, omega, box);
        integrator.getMoveManager().addMCMove(atomMove);




        getController().addActivity(new ActivityIntegrate(integrator));



        // make an IntegratorMC and write a move with a doTrial that samples the harmonic part of your system


    }

    public static void main(String[] args) {

        OctaneParams params = new OctaneParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            double temperature = 1;
            int nBeads = 10;
            boolean graphics = false;
            double omega = 1.0; // m*w^2
            double k4 = 1.0;
            long numSteps = 1000000;
        }

        double temperature = params.temperature;
        double omega = params.omega;;
        double k4 = params.k4;;
        int nBeads = params.nBeads;
        boolean graphics = params.graphics;
        long numSteps = params.numSteps;
        System.out.println(numSteps+" steps");

        final SimQuantumAO sim = new SimQuantumAO(Space1D.getInstance(), nBeads, temperature, omega, k4);
        sim.integrator.reset();

        MeterPrimPI meterPrimPI = new MeterPrimPI(sim.pmBonding, sim.pcP1, sim.betaN);





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
        long blockSize = numSteps/numBlocks;
        if (blockSize == 0) blockSize = 1;
        System.out.println("block size "+blockSize);
        AccumulatorAverageFixed accumulator = new AccumulatorAverageFixed(blockSize);
        DataPumpListener accumulatorPump = new DataPumpListener(meterPrimPI, accumulator);
        sim.integrator.getEventManager().addListener(accumulatorPump);

        final long startTime = System.currentTimeMillis();

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps));

        //MeterTargetTP.openFW("x"+numMolecules+".dat");
        //MeterTargetTP.closeFW();

        DataGroup data = (DataGroup)accumulator.getData();
        IData dataErr = data.getData(accumulator.ERROR.index);
        IData dataAvg = data.getData(accumulator.AVERAGE.index);
        IData dataCorrelation = data.getData(accumulator.BLOCK_CORRELATION.index);
        double avg = dataAvg.getValue(0);
        double err = dataErr.getValue(0);
        double cor = dataCorrelation.getValue(0);

        System.out.println("Energy_prim: " + avg + " +/- " + err + " cor: " + cor);

        long endTime = System.currentTimeMillis();
        System.out.println("time: " + (endTime - startTime)/1000.0);



    }

    public static class OctaneParams extends ParameterBase {
        public double temperature = 1;
        public int nBeads = 10;
        public boolean graphics = false;
        public double omega = 1.0; // m*w^2
        public double k4 = 1.0;
        public long numSteps = 1000000;
    }
}
