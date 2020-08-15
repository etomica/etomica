/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;


import etomica.action.activity.ActivityIntegrate2;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.AccumulatorAverageCovariance;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountSteps;
import etomica.data.IData;
import etomica.graphics.ColorScheme;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMD;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncatedForceShifted;
import etomica.potential.Potential2SoftSpherical;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;

import java.awt.*;
import java.util.Arrays;


public class SimLJHTTISuperSFMD extends Simulation {

    public final CoordinateDefinition coordinateDefinition;
    public IntegratorMD integrator;

    public Box box;
    public Boundary boundary;
    public int[] nCells;
    public Basis basis;
    public Primitive primitive;
    public PotentialMasterList potentialMaster;
    public Potential2SoftSpherical potential;
    public SpeciesSpheresMono species;

    public SimLJHTTISuperSFMD(Space _space, int numAtoms, double density, double temperature, double rc, int[] seeds) {
        super(_space);
        if (seeds != null) {
            setRandom(new RandomMersenneTwister(seeds));
        }
        species = new SpeciesSpheresMono(this, space);
        species.setIsDynamic(true);
        addSpecies(species);

        potentialMaster = new PotentialMasterList(this, space);
        potentialMaster.lrcMaster().setEnabled(false);

        double L = Math.pow(4.0 / density, 1.0 / 3.0);
        int n = (int) Math.round(Math.pow(numAtoms / 4, 1.0 / 3.0));
        boundary = new BoundaryRectangularPeriodic(space, n * L);
        box = this.makeBox(boundary);
        box.setNMolecules(species, numAtoms);

        integrator = new IntegratorVelocityVerlet(potentialMaster, getRandom(), 0.005, temperature, box);
        integrator.setIsothermal(true);

        primitive = new PrimitiveCubic(space, n * L);

        nCells = new int[]{n, n, n};
        Basis basisFCC = new BasisCubicFcc();
        basis = new BasisBigCell(space, basisFCC, nCells);

        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(new int[]{1, 1, 1});

        potential = new P2LennardJones(space, 1.0, 1.0);
        potential = new P2SoftSphericalTruncatedForceShifted(space, potential, rc);
        AtomType sphereType = species.getLeafType();
        potentialMaster.addPotential(potential, new AtomType[]{sphereType, sphereType});

        int cellRange = 2;
        potentialMaster.setRange(1.2 * rc);
        potentialMaster.setCellRange(cellRange);

        this.getController2().addActivity(new ActivityIntegrate2(integrator));

        integrator.getEventManager().addListener(potentialMaster.getNeighborManager(box));
    }

    /**
     * @param args filename containing simulation parameters
     */
    public static void main(String[] args) {
        //set up simulation parameters
        SimOverlapParam params = new SimOverlapParam();
        if (args.length == 0) {
            params.numAtoms = 500;
            params.numSteps = 10000;
            params.temperature = 1;
            params.density = 1;
            params.rc = 3;
            params.bpharm = 9;
            params.doD2 = true;
        } else {
            ParseArgs.doParseArgs(params, args);
        }
        double density = params.density;
        long numSteps = params.numSteps;
        final int numAtoms = params.numAtoms;
        double temperature = params.temperature;
        double rc = params.rc;
        double bpharm = params.bpharm;
        int[] seeds = params.randomSeeds;

        System.out.println("Running Lennard-Jones simulation");
        System.out.println(numAtoms + " atoms at density " + density + " and temperature " + temperature);
        System.out.println(numSteps + " steps");

        //instantiate simulation
        final SimLJHTTISuperSFMD sim = new SimLJHTTISuperSFMD(Space.getInstance(3), numAtoms, density, temperature, rc * Math.pow(density, -1.0 / 3.0), seeds);
        if (seeds == null) {
            seeds = ((RandomMersenneTwister) sim.getRandom()).getSeedArray();
        }
        System.out.println("Random seeds: " + Arrays.toString(seeds));
        if (false) {
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            simGraphic.setPaintInterval(sim.box, 1);
            ColorScheme colorScheme = new ColorScheme() {
                protected Color[] allColors;

                public Color getAtomColor(IAtom a) {
                    if (allColors == null) {
                        allColors = new Color[768];
                        for (int i = 0; i < 256; i++) {
                            allColors[i] = new Color(255 - i, i, 0);
                        }
                        for (int i = 0; i < 256; i++) {
                            allColors[i + 256] = new Color(0, 255 - i, i);
                        }
                        for (int i = 0; i < 256; i++) {
                            allColors[i + 512] = new Color(i, 0, 255 - i);
                        }
                    }
                    return allColors[(2 * a.getLeafIndex()) % 768];
                }
            };
            simGraphic.getDisplayBox(sim.box).setColorScheme(colorScheme);

            DisplayTextBox timer = new DisplayTextBox();
            DataSourceCountSteps counter = new DataSourceCountSteps(sim.integrator);
            DataPumpListener counterPump = new DataPumpListener(counter, timer, 100);
            sim.integrator.getEventManager().addListener(counterPump);
            simGraphic.getPanel().controlPanel.add(timer.graphic());

            simGraphic.makeAndDisplayFrame("LJ FCC");

            return;
        }

        //start simulation

        double L = Math.pow(numAtoms, 1.0 / 3.0);

        // meter needs lattice energy, so make it now
        sim.integrator.reset();
        MeterSolidDA meterSolid = new MeterSolidDA(sim.getSpace(), sim.potentialMaster, sim.coordinateDefinition, params.doD2);
        meterSolid.setTemperature(temperature);
        meterSolid.setPRes(temperature * bpharm);
        IData d = meterSolid.getData();
        double uLat = d.getValue(0);
        System.out.println("uLat: " + uLat);
        double pLat = d.getValue(1) - temperature * density;
        System.out.println("pLat: " + pLat);

        if (args.length == 0) {
            // quick initialization
            sim.initialize(numSteps / 10);
        } else {
            long nSteps = 50 + numAtoms * 3;
            sim.initialize(nSteps);
        }

        int numBlocks = 100;
        int interval = 10;
        long blockSize = numSteps / (numBlocks * interval);
        if (blockSize == 0) blockSize = 1;
        int o = 2;
        while (blockSize < numSteps / 5 && (numSteps != numBlocks * interval * blockSize)) {
            interval = 2 + (o % 2 == 0 ? (o / 2) : -(o / 2));
            if (interval < 1 || interval > numSteps / 5) {
                throw new RuntimeException("oops interval " + interval);
            }
            blockSize = numSteps / (numBlocks * interval);
            if (blockSize == 0) blockSize = 1;
            o++;
        }
        if (numSteps != numBlocks * interval * blockSize) {
            throw new RuntimeException("unable to find appropriate intervals");
        }
        System.out.println("block size " + blockSize + " interval " + interval);

        final AccumulatorAverageCovariance avgSolid = new AccumulatorAverageCovariance(blockSize);
        DataPumpListener pumpPU = new DataPumpListener(meterSolid, avgSolid, interval);
        sim.integrator.getEventManager().addListener(pumpPU);

        final long startTime = System.currentTimeMillis();

        sim.getController2().runActivityBlocking(new ActivityIntegrate2(sim.integrator), numSteps);
        long endTime = System.currentTimeMillis();
        System.out.println();

        IData avgRawData = avgSolid.getData(avgSolid.AVERAGE);
        IData errRawData = avgSolid.getData(avgSolid.ERROR);
        IData corRawData = avgSolid.getData(avgSolid.BLOCK_CORRELATION);

        double avgU = avgRawData.getValue(0);
        double errU = errRawData.getValue(0);
        double corU = corRawData.getValue(0);

        double avgP = avgRawData.getValue(1);
        double errP = errRawData.getValue(1);
        double corP = corRawData.getValue(1);

        double avgBUc = avgRawData.getValue(2);
        double errBUc = errRawData.getValue(2);
        double corBUc = corRawData.getValue(2);

        double avgZc = avgRawData.getValue(3);
        double errZc = errRawData.getValue(3);
        double corZc = corRawData.getValue(3);

        double avgDadv2 = avgRawData.getValue(4);
        double errDadv2 = errRawData.getValue(4);
        double corDadv2 = corRawData.getValue(4);


        System.out.print(String.format("Uraw:  % 21.15e  %10.4e  % 6.3f\n", avgU, errU, corU));
        System.out.print(String.format("Praw:  % 21.15e  %10.4e  % 6.3f\n", avgP, errP, corP));
        System.out.print(String.format("bUc:   % 21.15e  %10.4e  % 6.3f\n", avgBUc, errBUc, corBUc));
        System.out.print(String.format("Zc:    % 21.15e  %10.4e  % 6.3f\n", avgZc, errZc, corZc));
        System.out.print(String.format("P:     % 21.15e  %10.4e\n", avgZc * density * temperature + pLat + bpharm * temperature, errZc * density * temperature));

        if (params.doD2) {
            IData covData = avgSolid.getData(avgSolid.BLOCK_COVARIANCE);
            double y = avgU * numAtoms - uLat * numAtoms - 1.5 * temperature * (numAtoms - 1);
            double ey = errU * numAtoms;
            double avgCv = (avgRawData.getValue(5) - y * y) / (temperature * temperature);
            int n = avgRawData.getLength();
            double coru2u = covData.getValue(0 * n + 5) / Math.sqrt(covData.getValue(0 * n + 0) * covData.getValue(5 * n + 5));
            double errCv = Math.sqrt(errRawData.getValue(5) * errRawData.getValue(5) + 4 * y * y * ey * ey - 4 * y * ey * errRawData.getValue(5) * coru2u) / (temperature * temperature);
            double corCv = corRawData.getValue(5);

            y = avgBUc * numAtoms;
            ey = errBUc * numAtoms;
//            System.out.println("Cvcraw: "+avgRawData.getValue(6)/numAtoms+" "+errRawData.getValue(6)/numAtoms);
            double avgCvc = avgRawData.getValue(6) - y * y;
            coru2u = covData.getValue(2 * n + 6) / Math.sqrt(covData.getValue(2 * n + 2) * covData.getValue(6 * n + 6));
            double errCvc = Math.sqrt(errRawData.getValue(6) * errRawData.getValue(6) + 4 * y * y * ey * ey - 4 * y * ey * errRawData.getValue(6) * coru2u);
            double corCvc = corRawData.getValue(6);

            System.out.print(String.format("Cvraw: % 21.15e  %10.4e  % 6.3f\n", avgCv / numAtoms, errCv / numAtoms, corCv));
//            System.out.print(String.format("Cvc0:  % 21.15e  %10.4e  % 6.3f\n", avgRawData.getValue(6) / numAtoms, errCvc / numAtoms, corCvc));
            System.out.print(String.format("Cvc:   % 21.15e  %10.4e  % 6.3f\n", avgCvc / numAtoms, errCvc / numAtoms, corCvc));

            y = avgDadv2 * numAtoms;
            ey = errDadv2 * numAtoms;
            double avgCvc2 = avgRawData.getValue(7) - y * y;
            coru2u = covData.getValue(4 * n + 7) / Math.sqrt(covData.getValue(4 * n + 4) * covData.getValue(7 * n + 7));
            double errCvc2 = Math.sqrt(errRawData.getValue(7) * errRawData.getValue(7) + 4 * y * y * ey * ey - 4 * y * ey * errRawData.getValue(7) * coru2u);
            double corCvc2 = corRawData.getValue(7);

//            System.out.print(String.format("Cvcraw:% 21.15e  %10.4e  % 6.3f\n", avgRawData.getValue(7) / numAtoms, errRawData.getValue(7) / numAtoms, corCvc2));
            System.out.print(String.format("Cvc:   % 21.15e  %10.4e  % 6.3f\n", avgCvc2 / numAtoms, errCvc2 / numAtoms, corCvc2));

        }

        System.out.println();
        System.out.println("time: " + (endTime - startTime) / 1000.0);
    }

    public void initialize(long initSteps) {
        // equilibrate off the lattice to avoid anomalous contributions
        this.getController2().runActivityBlocking(new ActivityIntegrate2(this.integrator), initSteps);

    }

    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int numAtoms = 256;
        public double density = 1.28;
        public long numSteps = 100000;
        public double temperature = 0.1;
        public double rc = 2.5;
        public double bpharm = 0;
        public int[] randomSeeds = null;
        public boolean doD2 = false;
    }
}
