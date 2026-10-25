/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.normalmode;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.AccumulatorAverageCovariance;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountSteps;
import etomica.data.IData;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressure;
import etomica.data.types.DataGroup;
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
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;

import java.awt.*;

public class SimLJPropsCij extends Simulation {
    public final CoordinateDefinition coordinateDefinition;
    public IntegratorMD integrator;
    public Box box;
    public Boundary boundary;
    public int[] nCells;
    public Basis basis;
    public Primitive primitive;
    public PotentialMasterList potentialMaster;
    public SpeciesGeneral species;

    public SimLJPropsCij(Space _space, int N, double rho0, double T, double rc, double dt) {
        super(_space);

        setRandom(new RandomMersenneTwister(1)); // random seed

        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(species);
        double L = Math.pow(4.0/rho0, 1.0 / 3.0);
        int n = (int) Math.round(Math.pow(N / 4, 1.0 / 3.0));
        boundary = new BoundaryRectangularPeriodic(space, n * L);
        box = this.makeBox(boundary);
        box.setNMolecules(species, N);
        potentialMaster = new PotentialMasterList(getSpeciesManager(), box, 2, rc, BondingInfo.noBonding());
        potentialMaster.doAllTruncationCorrection = false;
        integrator = new IntegratorVelocityVerlet(potentialMaster, getRandom(), dt, T, box);
        integrator.setIsothermal(true);
        primitive = new PrimitiveCubic(space, n * L);
        nCells = new int[]{n, n, n};
        Basis basisFCC = new BasisCubicFcc();
        basis = new BasisBigCell(space, basisFCC, nCells);
        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(new int[]{1, 1, 1});

        IPotential2 potential = new P2LennardJones();
        potential = new P2SoftSphericalTruncated(potential, rc);
        AtomType sphereType = species.getLeafType();
        potentialMaster.setPairPotential(sphereType, sphereType, potential);
        potentialMaster.init();
        integrator.getEventManager().removeListener(potentialMaster);
        this.getController().addActivity(new ActivityIntegrate(integrator));
        ((P2SoftSphericalTruncated) potential).setTruncationRadius(0.6 * boundary.getBoxSize().getX(0));
    }

    public void initialize(long initSteps) {
        this.getController().runActivityBlocking(new ActivityIntegrate(this.integrator, initSteps));
    }

    /**
     * @param args filename containing simulation parameters
     */
    public static void main(String[] args) {
        final long startTime = System.currentTimeMillis();
        SimOverlapParam params = new SimOverlapParam();
        ParseArgs.doParseArgs(params, args);
        boolean isGraphic = params.isGraphic;;
        double rho0 = params.rho0;
        long steps = params.steps;
        long stepsEq = params.steps/10;
        long stepsShort = steps/10;

        double dt = params.dt;
        final int N = params.N;

        double T = params.T;
        double T2 = T * T;
        double rc = params.rc;
        double strain_x = params.strain_x;
        double strain_yz = params.strain_yz;
        double[] elasticParams = new double[2];
        elasticParams[0] = params.gV;
        elasticParams[1] = params.gVV;

        final SimLJPropsCij sim = new SimLJPropsCij(Space.getInstance(3), N, rho0, T, rc*Math.pow(rho0, -1.0 / 3.0), dt);
        double V = sim.box.getBoundary().volume();
        double rho = N / V;


        //Ulat
        sim.integrator.reset();//uLat
        MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        double uLat = meterPE.getDataAsScalar();

        //Plat
        MeterPressure meterP = new MeterPressure(sim.box, sim.potentialMaster);
        meterP.setTemperature(0);//ZERO means NO ideal gas component!
        double pLat = meterP.getDataAsScalar();
        meterP.setTemperature(T);


        System.out.println("N " + N + " V " + V + " L " + Math.pow(V, 1.0 / 3.0));
        System.out.println("rho " + rho + " T " + T);
        System.out.println("rc " + rc);
        System.out.println("steps " + steps + " stepsEq " + stepsEq + " stepsShort " + stepsShort);
        System.out.println("dt " + dt);

        System.out.println("\nHMA Elastic Parameters");
        System.out.println("gV  " + params.gV + " gVV  " + params.gVV);

        System.out.println("\nLattice properties");
        System.out.println("uLat  " + uLat);
        System.out.println("pLat  " + pLat);

        if (isGraphic) {
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
            simGraphic.makeAndDisplayFrame("MD FCC LJ");
            return;
        }

        boolean doD2 = true;
        MeterSolidHMALJ meterElastic = new MeterSolidHMALJ(sim.getSpace(), sim.potentialMaster, sim.coordinateDefinition, elasticParams, T, doD2);
        meterElastic.setTemperature(T);

        sim.initialize(stepsEq);

        int numBlocks = 100;
        int interval = 10;
        long blockSize = steps / (numBlocks * interval);
        if (blockSize == 0) blockSize = 1;
        int o = 2;
        while (blockSize < steps / 5 && (steps != numBlocks * interval * blockSize)) {
            interval = 2 + (o % 2 == 0 ? (o / 2) : -(o / 2));
            if (interval < 1 || interval > steps / 5) {
                throw new RuntimeException("oops interval " + interval);
            }
            blockSize = steps / (numBlocks * interval);
            if (blockSize == 0) blockSize = 1;
            o++;
        }
        if (steps != numBlocks * interval * blockSize) {
            throw new RuntimeException("unable to find appropriate intervals");
        }
        System.out.println("block size " + blockSize + "  interval " + interval);

        final AccumulatorAverageCovariance accumulatorElastic = new AccumulatorAverageCovariance(blockSize);

        DataPumpListener pumpElastic = new DataPumpListener(meterElastic, accumulatorElastic, interval);
        sim.integrator.getEventManager().addListener(pumpElastic);

        System.out.println("Short sim for Covariance  " + stepsShort + " steps");
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, stepsShort));

        IData dataElasticAvg0 = accumulatorElastic.getData(accumulatorElastic.AVERAGE);
        IData dataElasticErr0 = accumulatorElastic.getData(accumulatorElastic.ERROR);
        double uShift = dataElasticAvg0.getValue(1);
        double errUshift = dataElasticErr0.getValue(1);
        double pShift = dataElasticAvg0.getValue(3);
        double errPshift = dataElasticErr0.getValue(3);

        accumulatorElastic.reset();
        System.out.println("uShift  " + uShift + "  " + errUshift);
        System.out.println("pShift  " + pShift + "  " + errPshift);

        meterElastic.setShift(uShift, pShift);

        //Production
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps));

        IData dataElasticAvg = accumulatorElastic.getData(accumulatorElastic.AVERAGE);
        IData dataElasticErr = accumulatorElastic.getData(accumulatorElastic.ERROR);
        IData dataElasticCor = accumulatorElastic.getData(accumulatorElastic.BLOCK_CORRELATION);
        IData dataElasticCov = accumulatorElastic.getData(accumulatorElastic.COVARIANCE);

        /** First Derivatives*/
        //Bulk
        double Uc     = dataElasticAvg.getValue(0) + uShift;
        double errUc  = dataElasticErr.getValue(0);
        double corUc  = dataElasticCor.getValue(0);
        double Um      = dataElasticAvg.getValue(1) + uShift;
        double errUm   = dataElasticErr.getValue(1);
        double corUm   = dataElasticCor.getValue(1);
        double Pc     = dataElasticAvg.getValue(2) + pShift;
        double errPc  = dataElasticErr.getValue(2);
        double corPc  = dataElasticCor.getValue(2);
        double Pm      = dataElasticAvg.getValue(3) + pShift;
        double errPm   = dataElasticErr.getValue(3);
        double corPm   = dataElasticCor.getValue(3);

/** Second Derivatives*/
        int nd = dataElasticAvg.getLength();
        double varUc  = dataElasticCov.getValue(0 * (nd + 1)) * N * N;
        double varUm   = dataElasticCov.getValue(1 * (nd + 1)) * N * N;
        double varPc  = dataElasticCov.getValue(2 * (nd + 1));
        double varPm   = dataElasticCov.getValue(3 * (nd + 1));
        double covUcPc  = dataElasticCov.getValue(0 * nd + 2) * N;
        double covUmPm   = dataElasticCov.getValue(1 * nd + 3) * N;

//Bulk
        double Cvc = varUc / T2 / N;
        double Cvm  = dataElasticAvg.getValue(4)  +  varUm / T2 / N;

        double Bc  = dataElasticAvg.getValue(5) - V / T * varPc;
        double Bm   = dataElasticAvg.getValue(6) - V / T * varPm;

        double Gc = rho + covUcPc / T2;
        double Gm = dataElasticAvg.getValue(7) + covUmPm / T2;


//Print
        System.out.println("\n++++ First Derivatives ++++");
        System.out.println("Uc   " + Uc + " +/- " + errUc + "   cor " + corUc);
        System.out.println("Um    " + Um + " +/- " + errUm + "   cor " + corUm);
        System.out.println("Pc   " + Pc + " +/- " + errPc + "   cor " + corPc);
        System.out.println("Pm    " + Pm + " +/- " + errPm + "   cor " + corPm);

        System.out.println("\n ++++ Second Derivatives ++++");
        System.out.println("Cvc " + Cvc);
        System.out.println("Cvc " + Cvm);
        System.out.println("Bc  " + Bc);
        System.out.println("Bm  " + Bm);
        System.out.println("Gc  " + Gc);
        System.out.println("Gm  " + Gm);

        long endTime = System.currentTimeMillis();
        System.out.println("\ntime (min)  " + (endTime - startTime) / 1000.0/60.0);
    }


    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int N = 500;
        public double rho0 = 1.0;
        public long steps = 10000;
        public double dt = 0.002;
        public double T = 1;
        public double rc = 3.0;
        public boolean isGraphic = false;
        public double strain_x = 0.0;
        public double strain_yz = 0.0;
        public double gV = 3.094707617816900;
        public double gVV = -0.199514733065326;
    }
}