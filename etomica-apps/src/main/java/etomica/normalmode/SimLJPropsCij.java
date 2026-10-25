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
//        elasticParams[2] = params.gx1;
//        elasticParams[3] = params.gy1;
//        elasticParams[4] = params.gy4;
//        elasticParams[5] = params.gx11;
//        elasticParams[6] = params.gy11;
//        elasticParams[7] = params.gx44;
//        elasticParams[8] = params.gy44;
//        elasticParams[9] = params.gx12;
//        elasticParams[10] = params.gz12;


        final SimLJPropsCij sim = new SimLJPropsCij(Space.getInstance(3), N, rho0, T, rc*Math.pow(rho0, -1.0 / 3.0), dt);
        double V = sim.box.getBoundary().volume();
        double rho = N / V;

        System.out.println("N " + N + " V " + V + " L " + Math.pow(V, 1.0 / 3.0));
        System.out.println("rho " + rho + " T " + T);
        System.out.println("rc " + rc);
        System.out.println("steps " + steps + " stepsEq " + stepsEq + " stepsShort " + stepsShort);
        System.out.println("dt " + dt);

        System.out.println("\nHMA Elastic Parameters");
        System.out.println("gV  " + params.gV + " gVV  " + params.gVV);
//        System.out.println("gx1  " + params.gx1 + " gy1  " + params.gy1 + " gy4 " + params.gy4);
//        System.out.println("gx11  " + params.gx11 + " gy11  " + params.gy11 + " gx44 " + params.gx44);
//        System.out.println("gy44  " + params.gy44 + " gx12  " + params.gx12 + " gz12 " + params.gz12);


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


        sim.integrator.reset();
        double uLat = sim.integrator.getPotentialEnergy()/N;

        MeterPressure pMeterLat = new MeterPressure(sim.box, sim.integrator.getPotentialCompute());
        pMeterLat.setTemperature(0); // Lattice
        double pLat = pMeterLat.getDataAsScalar();

        System.out.println("\nLattice properties");
        System.out.println("uLat  " + uLat);
        System.out.println("pLat  " + pLat);

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

        IData avgRawData0 = accumulatorElastic.getData(accumulatorElastic.AVERAGE);
        IData errRawData0 = accumulatorElastic.getData(accumulatorElastic.ERROR);
        double uShift = avgRawData0.getValue(1);
        double errUshift = errRawData0.getValue(1);
        double pShift = avgRawData0.getValue(3);
        double errPshift = errRawData0.getValue(3);

        accumulatorElastic.reset();
        System.out.println("uShift  " + uShift + "  " + errUshift);
        System.out.println("pShift  " + pShift + "  " + errPshift);

        meterElastic.setShift(uShift, pShift);

        //Production
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps));
        System.out.println();

        IData avgRawData = accumulatorElastic.getData(accumulatorElastic.AVERAGE);
        IData errRawData = accumulatorElastic.getData(accumulatorElastic.ERROR);
        IData corRawData = accumulatorElastic.getData(accumulatorElastic.BLOCK_CORRELATION);
        IData covRawData = accumulatorElastic.getData(accumulatorElastic.COVARIANCE);

        /** First Derivatives*/
        //Bulk
        double Uc     = avgRawData.getValue(0) + uShift;
        double errUc  = errRawData.getValue(0);
        double corUc  = corRawData.getValue(0);
        double Um      = avgRawData.getValue(1) + uShift;
        double errUm   = errRawData.getValue(1);
        double corUm   = corRawData.getValue(1);
        double Pc     = avgRawData.getValue(2) + pShift;
        double errPc  = errRawData.getValue(2);
        double corPc  = corRawData.getValue(2);
        double Pm      = avgRawData.getValue(3) + pShift;
        double errPm   = errRawData.getValue(3);
        double corPm   = corRawData.getValue(3);
        //Normal stress
//        double p1Shift = - pShift;
//        double P1c     = avgRawData.getValue(4) + p1Shift;
//        double errP1c  = errRawData.getValue(4);
//        double corP1c  = corRawData.getValue(4);
//        double P1m      = avgRawData.getValue(5) + p1Shift;
//        double errP1m   = errRawData.getValue(5);
//        double corP1m   = corRawData.getValue(5);
//        double P2c     = avgRawData.getValue(6) + p1Shift;
//        double errP2c  = errRawData.getValue(6);
//        double corP2c  = corRawData.getValue(6);
//        double P2m      = avgRawData.getValue(7) + p1Shift;
//        double errP2m   = errRawData.getValue(7);
//        double corP2m   = corRawData.getValue(7);
//        double P3c     = avgRawData.getValue(8) + p1Shift;
//        double errP3c  = errRawData.getValue(8);
//        double corP3c  = corRawData.getValue(8);
//        double P3m      = avgRawData.getValue(9) + p1Shift;
//        double errP3m   = errRawData.getValue(9);
//        double corP3m   = corRawData.getValue(9);
        //Shear stress
//        double P4c     = avgRawData.getValue(10);
//        double errP4c  = errRawData.getValue(10);
//        double corP4c  = corRawData.getValue(10);
//        double P4m      = avgRawData.getValue(11);
//        double errP4m   = errRawData.getValue(11);
//        double corP4m   = corRawData.getValue(11);
//        double P5c     = avgRawData.getValue(12);
//        double errP5c  = errRawData.getValue(12);
//        double corP5c  = corRawData.getValue(12);
//        double P5m      = avgRawData.getValue(13);
//        double errP5m   = errRawData.getValue(13);
//        double corP5m   = corRawData.getValue(13);
//        double P6c     = avgRawData.getValue(14);
//        double errP6c  = errRawData.getValue(14);
//        double corP6c  = corRawData.getValue(14);
//        double P6m      = avgRawData.getValue(15);
//        double errP6m   = errRawData.getValue(15);
//        double corP6m   = corRawData.getValue(15);

/** Second Derivatives*/
        int nd = avgRawData.getLength();
        double varUc  = covRawData.getValue(0 * (nd + 1)) * N * N;
        double varUm   = covRawData.getValue(1 * (nd + 1)) * N * N;
        double varPc  = covRawData.getValue(2 * (nd + 1));
        double varPm   = covRawData.getValue(3 * (nd + 1));
        double covUcPc  = covRawData.getValue(0 * nd + 2) * N;
        double covUmPm   = covRawData.getValue(1 * nd + 3) * N;

//        double varP1c = covRawData.getValue(4 * (nd + 1));
//        double varP1m  = covRawData.getValue(5 * (nd + 1));
//        double varP2c = covRawData.getValue(6 * (nd + 1));
//        double varP2m  = covRawData.getValue(7 * (nd + 1));
//        double varP3c = covRawData.getValue(8 * (nd + 1));
//        double varP3m  = covRawData.getValue(9 * (nd + 1));
//        double varP4c = covRawData.getValue(10 * (nd + 1));
//        double varP4m  = covRawData.getValue(11 * (nd + 1));
//        double varP5c = covRawData.getValue(12 * (nd + 1));
//        double varP5m  = covRawData.getValue(13 * (nd + 1));
//        double varP6c = covRawData.getValue(14 * (nd + 1));
//        double varP6m  = covRawData.getValue(15 * (nd + 1));
//
//        double covP1P2c = covRawData.getValue(4 * nd + 6);
//        double covP1P2m  = covRawData.getValue(5 * nd + 7);
//        double covP1P3c = covRawData.getValue(4 * nd + 8);
//        double covP1P3m  = covRawData.getValue(5 * nd + 9);
//        double covP2P3c = covRawData.getValue(6 * nd + 8);
//        double covP2P3m  = covRawData.getValue(7 * nd + 9);

//        double covUP1c = covRawData.getValue(0 * nd + 4);
//        double covUP1m  = covRawData.getValue(1 * nd + 5);
//        double covUP2c = covRawData.getValue(0 * nd + 6);
//        double covUP2m  = covRawData.getValue(1 * nd + 7);
//        double covUP3c = covRawData.getValue(0 * nd + 8);
//        double covUP3m  = covRawData.getValue(1 * nd + 9);

//Bulk
        double Cvc = varUc / T2 / N;
        double Cvm  = avgRawData.getValue(4)  +  varUm / T2 / N;
        double Bc  = avgRawData.getValue(5) - V / T * varPc;
        double Bm   = avgRawData.getValue(6) - V / T * varPm;
        double Gc = rho + covUcPc / T2;
        double Gm = avgRawData.getValue(7) + covUmPm / T2;

        //Cij
        //normal
//        double C11c = avgRawData.getValue(19) - V / T * varP1c;
//        double C11m  = avgRawData.getValue(20) - V / T * varP1m;
//        double C22c = avgRawData.getValue(21) - V / T * varP2c;
//        double C22m  = avgRawData.getValue(22) - V / T * varP2m;
//        double C33c = avgRawData.getValue(23) - V / T * varP3c;
//        double C33m  = avgRawData.getValue(24) - V / T * varP3m;

//        double C12c = avgRawData.getValue(25) - V / T * covP1P2c;
//        double C12m  = avgRawData.getValue(26) - V / T * covP1P2m;
//        double C13c = avgRawData.getValue(27) - V / T * covP1P3c;
//        double C13m  = avgRawData.getValue(28) - V / T * covP1P3m;
//        double C23c = avgRawData.getValue(29) - V / T * covP2P3c;
//        double C23m  = avgRawData.getValue(30) - V / T * covP2P3m;
        //shear
//        double C44c = avgRawData.getValue(31) - V / T * varP4c;
//        double C44m  = avgRawData.getValue(32) - V / T * varP4m;
//        double C55c = avgRawData.getValue(33) - V / T * varP5c;
//        double C55m  = avgRawData.getValue(34) - V / T * varP5m;
//        double C66c = avgRawData.getValue(35) - V / T * varP6c;
//        double C66m  = avgRawData.getValue(36) - V / T * varP6m;

        //b_mn
//        double b11c = -rho + N / T2 * covUP1c;
//        double b11m  = avgRawData.getValue(38) + N / T2 * covUP1m;
//        double b22c = -rho + N/ T2 * covUP2c;
//        double b22m  = avgRawData.getValue(39) + N / T2 * covUP2m;
//        double b33c = -rho + N / T2 * covUP3c;
//        double b33m  = avgRawData.getValue(40) + N / T2 * covUP3m;

//Print
        System.out.println("\n ++++ First Derivatives ++++");
        System.out.print(String.format("Uc % 21.15f err  %10.4e  cor  % 6.3f\n", Uc, errUc, corUc));
        System.out.print(String.format("Um  % 21.15f err  %10.4e  cor  % 6.3f\n", Um, errUm, corUm));
        System.out.print(String.format("Pc % 21.15f err  %10.4e  cor  % 6.3f\n", Pc, errPc, corPc));
        System.out.print(String.format("Pm  % 21.15f err  %10.4e  cor  % 6.3f\n", Pm, errPm, corPm));
//        System.out.print(String.format(" P1c% 21.15f err  %10.4e  cor  % 6.3f\n", P1c, errP1c, corP1c));
//        System.out.print(String.format(" P1m % 21.15f err  %10.4e  cor  % 6.3f\n", P1m, errP1m, corP1m));
//        System.out.print(String.format(" P2c% 21.15f err  %10.4e  cor  % 6.3f\n", P2c, errP2c, corP2c));
//        System.out.print(String.format(" P2m % 21.15f err  %10.4e  cor  % 6.3f\n", P2m, errP2m, corP2m));
//        System.out.print(String.format(" P4c% 21.15f err  %10.4e  cor  % 6.3f\n", P4c, errP4c, corP4c));
//        System.out.print(String.format(" P4m % 21.15f err  %10.4e  cor  % 6.3f\n", P4m, errP4m, corP4m));


        System.out.println("\n ++++ Second Derivatives ++++");
        System.out.println("Cvc  " + Cvc);
        System.out.println("Cvc   " + Cvm);
        System.out.println("Bc   " + Bc);
        System.out.println("Bm    " + Bm);
        System.out.println("Gc " + Gc);
        System.out.println("Gm  " + Gm);

//        System.out.println("alphaVc " + gVc/Bc);
//        System.out.println("alphaVm  " + gVm/Bm);
//        System.out.println("Bsc  " + (Bc+T*V*gVc*gVc/Cvc));
//        System.out.println("Bsm   " + (Bm+T*V*gVm*gVm/Cvm));

//        System.out.println();
//        System.out.println("C11c     " + C11c);
//        System.out.println("C11c_avg " + 1.0/3.0*(C11c + C22c + C33c));
//        System.out.println("C11m      " + C11m);
//        System.out.println("C11m_avg  " + 1.0/3.0*(C11m + C22m + C33m));
//        System.out.println();
//        System.out.println("C12c     " + C12c);
//        System.out.println("C12c_avg " + 1.0/3.0*(C12c + C13c + C23c));
//        System.out.println("C12m      " + C12m);
//        System.out.println("C12m_avg  " + 1.0/3.0*(C12m + C13m + C23m));
//        System.out.println();
//        System.out.println("C44c     " + C44c);
//        System.out.println("C44c_avg " + 1.0/3.0*(C44c + C55c + C66c));
//        System.out.println("C44m      " + C44m);
//        System.out.println("C44m_avg  " + 1.0/3.0*(C44m + C55m + C66m));
//        System.out.println();
//        System.out.println("b11c     " + b11c);
//        System.out.println("b22c     " + b22c);
//        System.out.println("b11c_avg " + 1.0/3.0*(b11c + b22c + b33c));
//        System.out.println("b11m      " + b11m);
//        System.out.println("b22m      " + b22m);
//        System.out.println("b11m_avg  " + 1.0/3.0*(b11m + b22m + b33m));
//
        //Adiabatic Cij
//        double CijCorc = T*V/Cvc*gVc*gVc;
//        double CijCorm  = T*V/Cvm*gVc*gVc;
//        double C11sc = C11c + CijCorc;
//        double C11sm  = C11m  + CijCorm;
//        double C22sc = C22c + CijCorc;
//        double C22sm  = C22m  + CijCorm;
//        double C33sc = C33c + CijCorc;
//        double C33sm  = C33m  + CijCorm;
//        double C12sc = C12c + CijCorc;
//        double C12sm  = C12m  + CijCorm;
//        double C13sc = C13c + CijCorc;
//        double C13sm  = C13m  + CijCorm;
//        double C23sc = C23c + CijCorc;
//        double C23sm  = C23m  + CijCorm;

//        System.out.println("\n Adiabatic Cij");
//        System.out.println("C11sc     " + C11sc);
//        System.out.println("C11sc_avg " + 1.0/3.0*(C11sc + C22sc + C33sc));
//        System.out.println("C11sm      " + C11sm);
//        System.out.println("C11sm_avg  " + 1.0/3.0*(C11sm + C22sm + C33sm));
//        System.out.println();
//        System.out.println("C12sc     " + C12sc);
//        System.out.println("C12sc_avg " + 1.0/3.0*(C12sc + C13sc + C23sc));
//        System.out.println("C12sm      " + C12sm);
//        System.out.println("C12sm_avg  " + 1.0/3.0*(C12sm + C13sm + C23sm));
//        System.out.println();

        long endTime = System.currentTimeMillis();
        System.out.println("\ntime (min)  " + (endTime - startTime) / 1000.0/60.0);
    }


    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int N = 256;
        public double rho0 = 1.0;
        public long steps = 100000;
        public double dt = 0.002;
        public double T = 0.5;
        public double rc = 3.0;
        public boolean isGraphic = false;
        public double strain_x = 0.0;
        public double strain_yz = 0.0;

        //LJ: rc=3, rho=1.0
        public double gV = 3.094707617816900;
        public double gVV = -0.199514733065326;


//        public double gx1=3.577219390645709;
//        public double gy1=2.853451732227192;
//        public double gy4=2.232143500352495;
//        public double gx11=15.124623202977329;
//        public double gy11=4.801364959766065;
//        public double gx44=10.768197646340099;
//        public double gy44=3.233340182220370;
//        public double gx12=-4.113585532872547;
//        public double gz12=-0.392260472545948;
    }
}