/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.normalmode;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.chem.elements.Copper;
import etomica.data.AccumulatorAverageCovariance;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountSteps;
import etomica.data.IData;
import etomica.data.meter.MeterPressure;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorScheme;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.*;
import etomica.potential.compute.PotentialComputeEAM;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.units.*;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;

import java.awt.*;

/**
 * LREP Cu simulation with an FCC crystal, using HMA elastic method.
 */
public class SimLREPHMA extends Simulation {
    public final CoordinateDefinition coordinateDefinition;
    public IntegratorVelocityVerlet integrator;
    public SpeciesGeneral species;
    public Box box;
    public Boundary boundary;
    public int[] nCells;
    public Basis basis;
    public Primitive primitive;
    public PotentialComputeEAM potentialMaster;

    public SimLREPHMA(int numAtoms, double temperature, double density0, double dt) {
        super(Space3D.getInstance());
//        setRandom(new RandomMersenneTwister(2)); // set random seed for debugging
        species = SpeciesGeneral.monatomic(space, AtomType.element(Copper.INSTANCE), true);
        addSpecies(species);
        double L = Math.pow(4.0/density0, 1.0 / 3.0);
        int n = (int) Math.round(Math.pow(numAtoms / 4, 1.0 / 3.0));
        boundary = new BoundaryRectangularPeriodic(space, n * L);
        System.out.println(" L (A): " + boundary.getBoxSize().getX(0));
        box = this.makeBox(boundary);
        NeighborListManager nbrs = new NeighborListManager(this.getSpeciesManager(), box, 2, 7.8, BondingInfo.noBonding());
        nbrs.setDoDownNeighbors(true);

        potentialMaster = new PotentialComputeEAM(getSpeciesManager(), box, nbrs);
        potentialMaster.doAllTruncationCorrection = false;
        integrator = new IntegratorVelocityVerlet(potentialMaster, random, dt*1.0E-3, temperature, box);
        integrator.setIsothermal(true);

        box.setNMolecules(species, numAtoms);
        primitive = new PrimitiveCubic(space, n * L);
        nCells = new int[]{n, n, n};
        Basis basisFCC = new BasisCubicFcc();
        basis = new BasisBigCell(space, basisFCC, nCells);
        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(new int[]{1, 1, 1});

        AtomType leafType = species.getLeafType();
        P2LREPV p2 = new P2LREPV();
        potentialMaster.setPairPotential(leafType, leafType, p2);
        P2LREPPhi pRho = new P2LREPPhi();
        potentialMaster.setRhoPotential(leafType, pRho);
        EmbeddingSqrt f = new EmbeddingSqrt(1);
        potentialMaster.setEmbeddingPotential(leafType, f);

        potentialMaster.init();
        integrator.getEventManager().removeListener(nbrs);
        this.getController().addActivity(new ActivityIntegrate(integrator)); // for graphics
    }

    public static void main(String[] args) {
        final long startTime = System.currentTimeMillis();
        SimParams params = new SimParams();
        ParseArgs.doParseArgs(params, args);
        int numAtoms = params.numAtoms;
        double density0 = params.density0;
        double temperature = Kelvin.UNIT.toSim(params.temperatureK);
        boolean doGraphics = params.doGraphics;;
        boolean doD2 = params.doD2;
        long numSteps = params.numSteps;
        long numStepsEq = params.numSteps/10;

        double dt = params.dt;
        double[] elasticParams = new double[13];
        elasticParams[0] = params.gV;
        elasticParams[1] = params.gVV;
        elasticParams[2] = params.gx1;
        elasticParams[3] = params.gy1;
        elasticParams[4] = params.gy4;
        elasticParams[5] = params.gx11;
        elasticParams[6] = params.gy11;
        elasticParams[7] = params.gx44;
        elasticParams[8] = params.gy44;
        elasticParams[9] = params.gx12;
        elasticParams[10] = params.gz12;
        elasticParams[11] = params.dP;
        elasticParams[12] = params.dB;

        System.out.println(" gV: " + params.gV + " gVV: " + params.gVV);
        System.out.println(" gx1: " + params.gx1 + " gy1: " + params.gy1 + " gy4 " + params.gy4);
        System.out.println(" gx11: " + params.gx11 + " gy11: " + params.gy11 + " gx44 " + params.gx44);
        System.out.println(" gy44: " + params.gy44 + " gx12: " + params.gx12 + " gz12 " + params.gz12);
        System.out.println(" dP: " + params.dP + " dB: "+ params.dB);


        SimLREPHMA sim = new SimLREPHMA(numAtoms, temperature, density0, dt);

        if (doGraphics) {
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

        double volume = sim.box.getBoundary().volume();
        double density = numAtoms/volume;
        System.out.println(" numAtoms: " + numAtoms + "  volume: " + volume + "  density0: " + density0 + " density: " + density);
        System.out.println(" T (K): " + params.temperatureK);

        sim.integrator.reset();
        double uLat = 1e-3*Constants.AVOGADRO*Joule.UNIT.fromSim(sim.integrator.getPotentialEnergy()/numAtoms);
        System.out.println(" uLat (kJ/mol): " + uLat);

        MeterPressure pMeterLat = new MeterPressure(sim.box, sim.integrator.getPotentialCompute());
        pMeterLat.setTemperature(0); // Lattice
        double pLat = 1.0E-9*Pascal.UNIT.fromSim(pMeterLat.getDataAsScalar());
        System.out.println(" pLatPMeter (GPa): " + pLat);


        MeterSolidHMA meterElasticLat = new MeterSolidHMA(sim.getSpace(), sim.potentialMaster, sim.coordinateDefinition, elasticParams, temperature, doD2);
        sim.potentialMaster.computeAll(true, meterElasticLat);
        double virial = sim.potentialMaster.getLastVirial();
        double virial2 = sim.potentialMaster.getLastVirial2();
        double[] virialXYZ = sim.potentialMaster.getLastVirialXYZ();
        double[] virial2XYZ = sim.potentialMaster.getLastVirial2XYZ();
        double Plat  = -1e-9*Pascal.UNIT.fromSim(virial/3.0/volume);
        double P1lat = -1e-9*Pascal.UNIT.fromSim(virialXYZ[0]/volume);
        double P4lat = -1e-9*Pascal.UNIT.fromSim(virialXYZ[3]/volume);
        double Blat = 1e-9*Pascal.UNIT.fromSim(virial2/volume);
        double C11Lat = 1e-9*Pascal.UNIT.fromSim(virial2XYZ[0]/volume);
        double C12Lat = 1e-9*Pascal.UNIT.fromSim(virial2XYZ[3]/volume);
        double C44Lat = 1e-9*Pascal.UNIT.fromSim(virial2XYZ[6]/volume);
        System.out.println(" Plat (GPa): " + Plat + " P1lat (GPa): " + P1lat +" P4lat (GPa): " + P4lat);
        System.out.println( " Blat (GPa): " + Blat + " C11lat (GPa): " + C11Lat + " C12lat (GPa): " + C12Lat + " C44lat (GPa): " + C44Lat);



        //Initialization
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numStepsEq));
        System.out.println("**Done with "+ numStepsEq +" steps of equilibration ... ");

        int numBlocks = 100;
        int interval = 20;
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
        System.out.println(" equilibration: " + numStepsEq + "  production: " + numSteps + "  dt(fs): " + dt);
        System.out.println(" block size " + blockSize + "  interval " + interval);

        final AccumulatorAverageCovariance accumulatorElastic = new AccumulatorAverageCovariance(blockSize);

        MeterSolidHMA meterElastic = new MeterSolidHMA(sim.getSpace(), sim.potentialMaster, sim.coordinateDefinition, elasticParams, temperature, doD2);
        meterElastic.setTemperature(temperature);
        DataPumpListener pumpElastic = new DataPumpListener(meterElastic, accumulatorElastic, interval);
        sim.integrator.getEventManager().addListener(pumpElastic);

        double uShift = 0, pShift = 0;
        if(doD2){
            long numStepsShort = numSteps/10;
            System.out.println(" Short sim for Covariance: " + numStepsShort + " numSteps");
            sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numStepsShort));
            System.out.println("** Done with "+ numStepsShort +" steps of short run ... ");

            IData avgRawData0 = accumulatorElastic.getData(accumulatorElastic.AVERAGE);
            IData errRawData0 = accumulatorElastic.getData(accumulatorElastic.ERROR);
            uShift = avgRawData0.getValue(1); double errUshift = errRawData0.getValue(1);
            pShift = avgRawData0.getValue(3); double errPshift = errRawData0.getValue(3);
            uShift = 1e-3*Constants.AVOGADRO*Joule.UNIT.fromSim(uShift);
            errUshift = 1e-3*Constants.AVOGADRO*Joule.UNIT.fromSim(errUshift);
            pShift = 1e-9*Pascal.UNIT.fromSim(pShift);
            errPshift = 1e-9*Pascal.UNIT.fromSim(errPshift);
            System.out.println(" uShift (kJ/mol):  " + uShift + " +/- " + errUshift);
            System.out.println(" pShift (GPa):  " + pShift + " +/- " + errPshift);
            meterElastic.setShift(uShift, pShift);
            accumulatorElastic.reset();
        }


        //Production
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps));
        System.out.println();

        DataGroup data = (DataGroup)accumulatorElastic.getData();
        IData avgRawData = data.getData(accumulatorElastic.AVERAGE.index);
        IData errRawData = data.getData(accumulatorElastic.ERROR.index);
        IData corRawData = data.getData(accumulatorElastic.BLOCK_CORRELATION.index);
        IData covRawData = data.getData(accumulatorElastic.COVARIANCE.index);

        /** First Derivatives*/
        //Bulk
        double u_conv     = 1e-3*Constants.AVOGADRO*Joule.UNIT.fromSim(avgRawData.getValue(0) + uShift); //kJ/mol
        double errU_conv  = 1e-3*Constants.AVOGADRO*Joule.UNIT.fromSim(errRawData.getValue(0));
        double corU_conv  = corRawData.getValue(0);
        double u_hma      = 1e-3*Constants.AVOGADRO*Joule.UNIT.fromSim(avgRawData.getValue(1) + uShift);//GPa
        double errU_hma   = 1e-3*Constants.AVOGADRO*Joule.UNIT.fromSim(errRawData.getValue(1));
        double corU_hma   = corRawData.getValue(1);
        double P_conv     = 1e-9*Pascal.UNIT.fromSim(avgRawData.getValue(2) + pShift); //GPa
        double errP_conv  = 1e-9*Pascal.UNIT.fromSim(errRawData.getValue(2));
        double corP_conv  = corRawData.getValue(2);
        double P_hma      = 1e-9*Pascal.UNIT.fromSim(avgRawData.getValue(3) + pShift);
        double errP_hma   = 1e-9*Pascal.UNIT.fromSim(errRawData.getValue(3));
        double corP_hma   = corRawData.getValue(3);
        //Normal stress
        double p1Shift = - pShift;
        double P1_conv     = 1e-9*Pascal.UNIT.fromSim(avgRawData.getValue(4) + p1Shift);
        double errP1_conv  = 1e-9*Pascal.UNIT.fromSim(errRawData.getValue(4));
        double corP1_conv  = corRawData.getValue(4);
        double P1_hma      = 1e-9*Pascal.UNIT.fromSim(avgRawData.getValue(5) + p1Shift);
        double errP1_hma   = 1e-9*Pascal.UNIT.fromSim(errRawData.getValue(5));
        double corP1_hma   = corRawData.getValue(5);
        double P2_conv     = 1e-9*Pascal.UNIT.fromSim(avgRawData.getValue(6) + p1Shift);
        double errP2_conv  = 1e-9*Pascal.UNIT.fromSim(errRawData.getValue(6));
        double corP2_conv  = corRawData.getValue(6);
        double P2_hma      = 1e-9*Pascal.UNIT.fromSim(avgRawData.getValue(7) + p1Shift);
        double errP2_hma   = 1e-9*Pascal.UNIT.fromSim(errRawData.getValue(7));
        double corP2_hma   = corRawData.getValue(7);
        //Shear stress
        double P4_conv     = 1e-9*Pascal.UNIT.fromSim(avgRawData.getValue(10));
        double errP4_conv  = 1e-9*Pascal.UNIT.fromSim(errRawData.getValue(10));
        double corP4_conv  = corRawData.getValue(10);
        double P4_hma      = 1e-9*Pascal.UNIT.fromSim(avgRawData.getValue(11));
        double errP4_hma   = 1e-9*Pascal.UNIT.fromSim(errRawData.getValue(11));
        double corP4_hma   = corRawData.getValue(11);

        System.out.println("\n ++++ First Derivatives ++++");
        System.out.print(String.format(" u_conv (kJ/mol): % 21.15f err: %10.4e  cor: % 6.3f\n", u_conv, errU_conv, corU_conv));
        System.out.print(String.format(" u_hma (kJ/mol):  % 21.15f err: %10.4e  cor: % 6.3f\n", u_hma, errU_hma, corU_hma));
        System.out.print(String.format(" P_conv (GPa):    % 21.15f err: %10.4e  cor: % 6.3f\n", P_conv, errP_conv, corP_conv));
        System.out.print(String.format(" P_hma (GPa):     % 21.15f err: %10.4e  cor: % 6.3f\n", P_hma, errP_hma, corP_hma));
        System.out.print(String.format(" P1_conv (GPa):   % 21.15f err: %10.4e  cor: % 6.3f\n", P1_conv, errP1_conv, corP1_conv));
        System.out.print(String.format(" P1_hma (GPa):    % 21.15f err: %10.4e  cor: % 6.3f\n", P1_hma, errP1_hma, corP1_hma));
        System.out.print(String.format(" P2_conv (GPa):   % 21.15f err: %10.4e  cor: % 6.3f\n", P2_conv, errP2_conv, corP2_conv));
        System.out.print(String.format(" P2_hma (GPa):    % 21.15f err: %10.4e  cor: % 6.3f\n", P2_hma, errP2_hma, corP2_hma));
        System.out.print(String.format(" P4_conv (GPa):   % 21.15f err: %10.4e  cor: % 6.3f\n", P4_conv, errP4_conv, corP4_conv));
        System.out.print(String.format(" P4_hma (GPa):    % 21.15f err: %10.4e  cor: % 6.3f\n", P4_hma, errP4_hma, corP4_hma));


        if(doD2){
            /** Second Derivatives*/
            int nd = avgRawData.getLength();//41
            double varU_conv  = covRawData.getValue(0 * (nd + 1));
            double varU_hma   = covRawData.getValue(1 * (nd + 1));
            double varP_conv  = covRawData.getValue(2 * (nd + 1));
            double varP_hma   = covRawData.getValue(3 * (nd + 1));
            double varP1_conv = covRawData.getValue(4 * (nd + 1));
            double varP1_hma  = covRawData.getValue(5 * (nd + 1));
            double varP2_conv = covRawData.getValue(6 * (nd + 1));
            double varP2_hma  = covRawData.getValue(7 * (nd + 1));
            double varP3_conv = covRawData.getValue(8 * (nd + 1));
            double varP3_hma  = covRawData.getValue(9 * (nd + 1));
            double varP4_conv = covRawData.getValue(10 * (nd + 1));
            double varP4_hma  = covRawData.getValue(11 * (nd + 1));
            double varP5_conv = covRawData.getValue(12 * (nd + 1));
            double varP5_hma  = covRawData.getValue(13 * (nd + 1));
            double varP6_conv = covRawData.getValue(14 * (nd + 1));
            double varP6_hma  = covRawData.getValue(15 * (nd + 1));

            double covP1P2_conv = covRawData.getValue(4 * nd + 6);
            double covP1P2_hma  = covRawData.getValue(5 * nd + 7);
            double covP1P3_conv = covRawData.getValue(4 * nd + 8);
            double covP1P3_hma  = covRawData.getValue(5 * nd + 9);
            double covP2P3_conv = covRawData.getValue(6 * nd + 8);
            double covP2P3_hma  = covRawData.getValue(7 * nd + 9);
//
            double covUP_conv  = covRawData.getValue(0 * nd + 2);
            double covUP_hma   = covRawData.getValue(1 * nd + 3);
            double covUP1_conv = covRawData.getValue(0 * nd + 4);
            double covUP1_hma  = covRawData.getValue(1 * nd + 5);
            double covUP2_conv = covRawData.getValue(0 * nd + 6);
            double covUP2_hma  = covRawData.getValue(1 * nd + 7);
            double covUP3_conv = covRawData.getValue(0 * nd + 8);
            double covUP3_hma  = covRawData.getValue(1 * nd + 9);

//Bulk
            double Cv_conv = numAtoms*numAtoms*varU_conv/temperature/temperature + 3.0/2.0*(numAtoms-1);
            double Cv_hma  = avgRawData.getValue(16) + numAtoms*numAtoms*varU_hma/temperature / temperature + 3.0*(numAtoms-1.0);
            double B_conv  = avgRawData.getValue(17) - volume/ temperature * varP_conv;
            double B_hma   = avgRawData.getValue(18) - volume / temperature * varP_hma;
//Cij
            //normal
            double C11_conv = avgRawData.getValue(19) - volume / temperature * varP1_conv;
            double C11_hma  = avgRawData.getValue(20) - volume / temperature * varP1_hma;
            double C22_conv = avgRawData.getValue(21) - volume / temperature * varP2_conv;
            double C22_hma  = avgRawData.getValue(22) - volume / temperature * varP2_hma;
            double C33_conv = avgRawData.getValue(23) - volume / temperature * varP3_conv;
            double C33_hma  = avgRawData.getValue(24) - volume / temperature * varP3_hma;
//
            double C12_conv = avgRawData.getValue(25) - volume / temperature * covP1P2_conv;
            double C12_hma  = avgRawData.getValue(26) - volume / temperature * covP1P2_hma;
            double C13_conv = avgRawData.getValue(27) - volume / temperature * covP1P3_conv;
            double C13_hma  = avgRawData.getValue(28) - volume / temperature * covP1P3_hma;
            double C23_conv = avgRawData.getValue(29) - volume / temperature * covP2P3_conv;
            double C23_hma  = avgRawData.getValue(30) - volume / temperature * covP2P3_hma;
            //shear
            double C44_conv = avgRawData.getValue(31) - volume / temperature * varP4_conv;
            double C44_hma  = avgRawData.getValue(32) - volume / temperature * varP4_hma;
            double C55_conv = avgRawData.getValue(33) - volume / temperature * varP5_conv;
            double C55_hma  = avgRawData.getValue(34) - volume / temperature * varP5_hma;
            double C66_conv = avgRawData.getValue(35) - volume / temperature * varP6_conv;
            double C66_hma  = avgRawData.getValue(36) - volume / temperature * varP6_hma;
//          //gV
            double gV_conv  = density + numAtoms / temperature / temperature * covUP_conv;
            double gV_hma   = avgRawData.getValue(37) + numAtoms / temperature / temperature * covUP_hma;
//            //b_mn
            double b11_conv = -density + numAtoms / temperature / temperature * covUP1_conv;
            double b11_hma  = avgRawData.getValue(38) + numAtoms / temperature / temperature * covUP1_hma;
            double b22_conv = -density + numAtoms/ temperature / temperature * covUP2_conv;
            double b22_hma  = avgRawData.getValue(39) + numAtoms / temperature / temperature * covUP2_hma;
            double b33_conv = -density + numAtoms / temperature / temperature * covUP3_conv;
            double b33_hma  = avgRawData.getValue(40) + numAtoms / temperature / temperature * covUP3_hma;
//
//Print

            System.out.println("\n ++++ Second Derivatives ++++");
            UnitRatio cvUnit = new UnitRatio(Joule.UNIT, Kelvin.UNIT);
            System.out.println(" Cv_conv (J/mol-K):   " + Constants.AVOGADRO*cvUnit.fromSim(Cv_conv/numAtoms)); // J/mol-K
            System.out.println(" Cv_hma  (J/mol-K):   " + Constants.AVOGADRO*cvUnit.fromSim(Cv_hma/numAtoms)); // J/mol-K
            System.out.println(" B_conv (GPa):        " + 1e-9*Pascal.UNIT.fromSim(B_conv));
            System.out.println(" B_hma  (GPa):        " + 1e-9*Pascal.UNIT.fromSim(B_hma));
            UnitRatio gvUnit = new UnitRatio(Pascal.UNIT, Kelvin.UNIT);
            System.out.println(" gV_conv (GPa/K):     " + 1e-9*gvUnit.fromSim(gV_conv));//GPa/K
            System.out.println(" gV_hma  (GPa/K):     " + 1e-9*gvUnit.fromSim(gV_hma));
            double alphaV_conv = gV_conv/B_conv;
            System.out.println(" alphaV_conv (K^-1):  " + 1.0/Kelvin.UNIT.fromSim(1.0/alphaV_conv));
            double alphaV_hma = gV_hma/B_hma;
            System.out.println(" alphaV_hma  (K^-1):  " + 1.0/Kelvin.UNIT.fromSim(1.0/alphaV_hma));


            System.out.println(" Bs_conv (GPa):       " + 1e-9*Pascal.UNIT.fromSim((B_conv + temperature * volume * gV_conv * gV_conv / Cv_conv)));
            System.out.println(" Bs_hma  (GPa):       " + 1e-9*Pascal.UNIT.fromSim((B_hma+temperature*volume*gV_hma*gV_hma/Cv_hma)));
//
            System.out.println();
            System.out.println(" C11_conv     (GPa):  "     + 1e-9*Pascal.UNIT.fromSim(C11_conv));
            System.out.println(" C11_conv_avg (GPa):  " + 1e-9*Pascal.UNIT.fromSim(1.0/3.0*(C11_conv + C22_conv + C33_conv)));

            System.out.println(" C11_hma     (GPa):   "      + 1e-9*Pascal.UNIT.fromSim(C11_hma));
            System.out.println(" C11_hma_avg (GPa):   "  + 1e-9*Pascal.UNIT.fromSim(1.0/3.0*(C11_hma + C22_hma + C33_hma)));

            System.out.println();
            System.out.println(" C12_conv     (GPa):  "     + 1e-9*Pascal.UNIT.fromSim(C12_conv));
            System.out.println(" C12_conv_avg (GPa):  " + 1e-9*Pascal.UNIT.fromSim(1.0/3.0*(C12_conv + C13_conv + C23_conv)));

            System.out.println(" C12_hma     (GPa):   "      + 1e-9*Pascal.UNIT.fromSim(C12_hma));
            System.out.println(" C12_hma_avg (GPa):   "  + 1e-9*Pascal.UNIT.fromSim(1.0/3.0*(C12_hma + C13_hma + C23_hma)));

            System.out.println();
            System.out.println(" C44_conv     (GPa):  "     + 1e-9*Pascal.UNIT.fromSim(C44_conv));
            System.out.println(" C44_conv_avg (GPa):  " + 1e-9*Pascal.UNIT.fromSim(1.0/3.0*(C44_conv + C55_conv + C66_conv)));

            System.out.println(" C44_hma     (GPa):   "      + 1e-9*Pascal.UNIT.fromSim(C44_hma));
            System.out.println(" C44_hma_avg (GPa):   "  + 1e-9*Pascal.UNIT.fromSim(1.0/3.0*(C44_hma + C55_hma + C66_hma)));
            System.out.println();

            System.out.println(" gV_conv      (GPa/K):  " + 1e-9*gvUnit.fromSim(gV_conv));//GPa/K
            System.out.println(" b11_conv     (GPa/K):  " + 1e-9*gvUnit.fromSim(b11_conv));
            System.out.println(" b11_conv_avg (GPa/K):  " + 1e-9*gvUnit.fromSim(1.0/3.0*(b11_conv + b22_conv + b33_conv)));
            System.out.println(" b11_hma      (GPa/K):  " + 1e-9*gvUnit.fromSim(b11_hma));
            System.out.println(" b11_hma_avg  (GPa/K):  " + 1e-9*gvUnit.fromSim(1.0/3.0*(b11_hma + b22_hma + b33_hma)));
//
            //Adiabatic Cij
            double CijCor_conv = temperature*volume/Cv_conv*gV_conv*gV_conv;
            double CijCor_hma  = temperature*volume/Cv_hma*gV_conv*gV_conv;
            double C11s_conv = C11_conv + CijCor_conv;
            double C11s_hma  = C11_hma  + CijCor_hma;
            double C22s_conv = C22_conv + CijCor_conv;
            double C22s_hma  = C22_hma  + CijCor_hma;
            double C33s_conv = C33_conv + CijCor_conv;
            double C33s_hma  = C33_hma  + CijCor_hma;
            double C12s_conv = C12_conv + CijCor_conv;
            double C12s_hma  = C12_hma  + CijCor_hma;
            double C13s_conv = C13_conv + CijCor_conv;
            double C13s_hma  = C13_hma  + CijCor_hma;
            double C23s_conv = C23_conv + CijCor_conv;
            double C23s_hma  = C23_hma  + CijCor_hma;

            System.out.println("\n Adiabatic Cij");
            System.out.println(" C11s_conv (GPa):     " + 1e-9*Pascal.UNIT.fromSim(C11s_conv));
            System.out.println(" C11s_conv_avg (GPa): " + 1e-9*Pascal.UNIT.fromSim(1.0/3.0*(C11s_conv + C22s_conv + C33s_conv)));
            System.out.println(" C11s_hma (GPa):      " + 1e-9*Pascal.UNIT.fromSim(C11s_hma));
            System.out.println(" C11s_hma_avg (GPa):  " + 1e-9*Pascal.UNIT.fromSim(1.0/3.0*(C11s_hma + C22s_hma + C33s_hma)));
            System.out.println();
            System.out.println(" C12s_conv (GPa):     " + 1e-9*Pascal.UNIT.fromSim(C12s_conv));
            System.out.println(" C12s_conv_avg (GPa): " + 1e-9*Pascal.UNIT.fromSim(1.0/3.0*(C12s_conv + C13s_conv + C23s_conv)));
            System.out.println(" C12s_hma (GPa):      " + 1e-9*Pascal.UNIT.fromSim(C12s_hma));
            System.out.println(" C12s_hma_avg (GPa):  " + 1e-9*Pascal.UNIT.fromSim(1.0/3.0*(C12s_hma + C13s_hma + C23s_hma)));
            System.out.println();
        }

        long endTime = System.currentTimeMillis();
        System.out.println("time (min): " + (endTime - startTime) / 1000.0/60.0);



    }

    public static class SimParams extends ParameterBase {
        public boolean doGraphics = false;
        public boolean doD2 = !true;

        public double dt = 1.0;
        public int numAtoms = 500;
        public long numSteps = 2000;
        public double temperatureK = 300;
        public double density0 = 0.08502338387498792;
        public double dP=1.63878;
        public double dB=-7.8;
        //EAM: DB - rho=0.08502338387498792
        public double gV=1.404500637647537;
        public double gVV=3.131965483325428;
        public double gx1=1.170601829131321;
        public double gy1=1.588168559881036;
        public double gx11=22.383238919513573;
        public double gy11=12.238653033941443;
        public double gy4=1.335764638983127;
        public double gx44=9.933754280102567;
        public double gy44=7.901601526408648;
        public double gx12=-0.721140067130724;
        public double gz12=-4.342584087701786;

        //EAM: FB - rho=0.08502338387498792
//        public double gV=1.665191135820476;
//        public double gVV=1.744027161038234;
//        public double gx1=2.218272165051940;
//        public double gy1=1.514827142811007;
//        public double gx11=10.263285896692766;
//        public double gy11=-2.828872977179797;
//        public double gy4=1.514827142815689;
//        public double gx44=2.538760288462838;
//        public double gy44-0.502663712521145;
//        public double gx12=1.628531514391176;
//        public double gz12=7.128162834496251;

    }
}