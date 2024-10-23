/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.config.ConformationLinear;
import etomica.data.*;
import etomica.data.types.DataFunction;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorScheme;
import etomica.graphics.DisplayBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorLangevin;
import etomica.integrator.IntegratorMD;
import etomica.lattice.LatticeCubicFcc;
import etomica.potential.*;
import etomica.potential.compute.PotentialCompute;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputePair;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * PIMD using Langevin Thermostat using BAOAB algorithm.
 *
 * @author Andrew Schultz
 */

public class LJPIMD2 extends Simulation {

    public final PotentialComputePair potentialMaster;
    public final IntegratorMD integrator;
    public final Box box;
    public final MCMoveHOReal2 moveStageSimple, moveStageEC;
    public double beta;
    public int dim;

    public LJPIMD2(Space space, MoveChoice coordType, double mass, double timeStep, double gammaLangevin, int nBeads, int numAtoms, double temperature, double density, double rc, double omega2, double hbar) {
        super(Space3D.getInstance());
        SpeciesGeneral species = new SpeciesBuilder(space)
                .setDynamic(true)
                .addCount(AtomType.simple("A", mass / nBeads), nBeads)
                .withConformation(new ConformationLinear(space, 0))
                .build();
        addSpecies(species);
        this.dim = space.D();
        box = new Box(space);
        addBox(box);
        NeighborListManagerPI neighborManager = new NeighborListManagerPI(getSpeciesManager(), box, 2, rc, BondingInfo.noBonding());
        neighborManager.setAutoUpdateNeighbors(false);
        potentialMaster = new PotentialComputePair(getSpeciesManager(), box, neighborManager);
        beta = 1.0/temperature;

        box.setNMolecules(species, numAtoms);
        BoxInflate inflater = new BoxInflate(box, space, density);
        inflater.actionPerformed();
        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        config.initializeCoordinates(box);

        P2LennardJones p2lj = new P2LennardJones(1.0, 1.0 / nBeads);
        P2SoftSphericalTruncated p2 = new P2SoftSphericalTruncated(p2lj, rc);
        AtomType atomType = species.getLeafType();
        potentialMaster.setPairPotential(atomType, atomType, p2);
        potentialMaster.doAllTruncationCorrection = false;

        PotentialComputeAggregate.localStorageDefault = true;

        if (coordType == LJPIMD2.MoveChoice.Real) {
            integrator = new IntegratorLangevin(potentialMaster, random, timeStep, temperature, box, gammaLangevin);
        } else if (coordType == LJPIMD2.MoveChoice.NM) {
            MCMoveHO move = new MCMoveHO(space, potentialMaster, random, temperature, 0, box, hbar);
            integrator = new IntegratorLangevinPINM(potentialMaster, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
        } else if (coordType == LJPIMD2.MoveChoice.NMEC) {
            MCMoveHO move = new MCMoveHO(space, potentialMaster, random, temperature, omega2, box, hbar);
            integrator = new IntegratorLangevinPINM(potentialMaster, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
        } else if (coordType == LJPIMD2.MoveChoice.Stage) {
            MCMoveHOReal2 move = new MCMoveHOReal2(space, potentialMaster, random, temperature, 0, box, hbar);
            integrator = new IntegratorLangevinPI(potentialMaster, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
        } else { //StageEC -- default
            MCMoveHOReal2 move = new MCMoveHOReal2(space, potentialMaster, random, temperature, omega2, box, hbar);
            integrator = new IntegratorLangevinPI(potentialMaster, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
        }

        // find neighbors now.  Don't hook up NeighborListManager (neighbors won't change)
        potentialMaster.init();
        // and then never update neighbor lists
        integrator.getEventManager().removeListener(neighborManager);
        p2.setTruncationRadius(2*rc);

        moveStageSimple = new MCMoveHOReal2(space, potentialMaster, random, temperature, 0, box, hbar);
        moveStageEC = new MCMoveHOReal2(space, potentialMaster, random, temperature, omega2, box, hbar);
        integrator.setThermostatNoDrift(false);
        integrator.setIsothermal(true);
    }

    public static void main(String[] args) {
        long t1 = System.currentTimeMillis();
        LJPIMD2.SimParams params = new LJPIMD2.SimParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.steps = 10000;
            params.hbar = 1;
            params.temperature = 0.5;
            params.numAtoms = 108;
            params.rc = 2.5;
            params.isGraphic = false;
            params.coordType = LJPIMD2.MoveChoice.StageEC;
            params.timeStep = 0.01;
            params.nBeads = 1;
        }
        double facTimestep = params.facTimestep;
        Space space = Space.getInstance(params.D);
        int nShifts = params.nShifts;
        double mass = params.mass;
        double hbar = params.hbar;
        int numAtoms = params.numAtoms;
        int nBeads = params.nBeads;
        double temperature = params.temperature;
        double density = params.density;
        double rc = params.rc;
        double k2 = params.k2;
        double omega2 = params.k2 / mass;
        double gammaLangevin = params.gammaLangevin;
        double timeStep = params.timeStep;
        boolean isGraphic = params.isGraphic;
        LJPIMD2.MoveChoice coordType = params.coordType;
        long steps = params.steps;
        long stepsEq = steps/10;

        double omega = Math.sqrt(omega2);
        double x = 1/temperature*hbar*omega;
        if (nBeads == -1){
            nBeads = (int) (20*x);
        }
        double omegaN = nBeads*temperature/hbar;

        LJPIMD2 sim = new LJPIMD2(space, coordType, mass, timeStep, gammaLangevin, nBeads, numAtoms, temperature, density, rc, omega2, hbar);
        sim.integrator.reset();

        System.out.println(" LJ PIMD-"+coordType);
        System.out.println(" mass: " + mass);
        System.out.println(" T: " + temperature);
        System.out.println(" hbar: " + hbar);
        System.out.println(" w: " + omega);
        System.out.println(" wn: " + omegaN  + " , w/sqrt(n): " + omega/Math.sqrt(nBeads));
        System.out.println(" x: beta*hbar*w = " + hbar*omega/temperature);
        System.out.println(" nBeads: " + nBeads);
        System.out.println(" nShifts: "+ nShifts);
        System.out.println(" steps: " +  steps + " stepsEq: " + stepsEq);
        System.out.println(" timestep: " + timeStep);
        System.out.println(" facTimestep: " + facTimestep);
        System.out.println(" k2: " + k2);
        System.out.println(" gammaLangevin: " + gammaLangevin);
        System.out.println(" N: " + numAtoms);
        System.out.println(" density: " + density);
        System.out.println(" rc: " + rc);


        MeterPIVir2 meterVir2 = new MeterPIVir2(sim.potentialMaster, temperature, sim.box);


        int interval = 5;
        int blocks = 100;
        long blockSize = steps / (interval * blocks);
        if (blockSize == 0) blockSize = 1;
        System.out.println(" numBlocks: " + blocks + " blocksize: " + blockSize + " interval: " + interval);

        if (isGraphic) {
            sim.getController().addActivity(new ActivityIntegrate(sim.integrator));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.GRAPHIC_ONLY);
            simGraphic.setPaintInterval(sim.box, 1);
            int finalNBeads = nBeads;

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
                    return allColors[(768 * a.getIndex() / (finalNBeads))];
                }
            };

            simGraphic.getDisplayBox(sim.box).setColorScheme(colorScheme);
            ((DiameterHashByType) ((DisplayBox) simGraphic.displayList().getFirst()).getDiameterHash()).setDiameter(sim.species().getAtomType(0), 1);
            simGraphic.makeAndDisplayFrame("PIMD-"+coordType);

            return;
        }

        System.out.flush();

        sim.potentialMaster.computeAll(true);
        double uLat = sim.potentialMaster.getLastEnergy()/numAtoms;
        double volume = sim.box.getBoundary().volume();
        double pLat = -sim.potentialMaster.getLastVirial()/3.0/volume;
        System.out.println(" uLat: " + uLat);
        System.out.println(" pLat: " + pLat);

        // equilibration
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, stepsEq));
        System.out.println(" equilibration finished");

        AccumulatorAverageCovariance accumulatorVir2 = new AccumulatorAverageCovariance(blockSize);
        DataPumpListener accumulatorPumpVir2 = new DataPumpListener(meterVir2, accumulatorVir2, interval);
        sim.integrator.getEventManager().addListener(accumulatorPumpVir2);

        //Short run
        double EnShift = 0, errEnShift = 0;
        long numStepsShort = steps/10;
        System.out.println(" Short sim for Covariance: " + numStepsShort + " numSteps");
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numStepsShort));
        System.out.println(" Done with "+ numStepsShort +" steps of short run");

        double beta = 1/temperature;
        DataGroup dataVir2 = (DataGroup) accumulatorVir2.getData();
        IData dataAvgVir2 = dataVir2.getData(accumulatorVir2.AVERAGE.index);
        IData dataErrVir2 = dataVir2.getData(accumulatorVir2.ERROR.index);
        IData dataCovVir2 = dataVir2.getData(accumulatorVir2.COVARIANCE.index);

        double xHx = dataAvgVir2.getValue(1);
        double sdot = -1.0/beta/(1.0+beta*xHx/(sim.dim*numAtoms));
        System.out.println(" sdot*beta: " + sdot*beta);
        meterVir2.setSdot(sdot);
        accumulatorVir2.reset();


        // Run
        sim.integrator.resetStepCount();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps));

        double FX2 = dataAvgVir2.getValue(2);
        double varFX = dataCovVir2.getValue(0);
        double uconv = dataAvgVir2.getValue(3);
        double uhma = dataAvgVir2.getValue(4);
        double uhma2 = dataAvgVir2.getValue(5);
        double errUconv = dataErrVir2.getValue(3);
        double errUhma = dataErrVir2.getValue(4);
        double errUhma2 = dataErrVir2.getValue(5);

        System.out.println();
        System.out.println(" uConv: " + uconv + " " + errUconv);
        System.out.println(" uHMA: " + uhma + " " + errUhma);
        System.out.println(" uHMA2: " + uhma2 + " " + errUhma2);


        long t2 = System.currentTimeMillis();
        System.out.println("\n time: (min) " + (t2 - t1) * 0.001 / 60.0);
    }

    public static void writeDataToFile(IDataSource meter, IData errData, String filename) throws IOException {
        IData data;
        IData xData;
        if (meter instanceof AccumulatorAverage) {
            AccumulatorAverage acc = (AccumulatorAverage) meter;
            data = acc.getData(acc.AVERAGE);
            xData = ((DataFunction.DataInfoFunction) ((DataGroup.DataInfoGroup) acc.getDataInfo()).getSubDataInfo(acc.AVERAGE.index)).getXDataSource().getIndependentData(0);
        } else {
            data = meter.getData();
            xData = ((DataFunction.DataInfoFunction) meter.getDataInfo()).getXDataSource().getIndependentData(0);
        }
        boolean allNaN = true;
        for (int i = 0; i < xData.getLength(); i++) {
            if (!Double.isNaN(data.getValue(i))) allNaN = false;
        }
        if (allNaN) return;
        FileWriter fw = new FileWriter(filename);
        for (int i = 0; i < xData.getLength(); i++) {
            double y = data.getValue(i);
            if (Double.isNaN(y)) continue;
            if (errData == null) {
                fw.write(xData.getValue(i) + " " + y + "\n");
            } else {
                fw.write(xData.getValue(i) + " " + y + " " + errData.getValue(i) + "\n");
            }
        }
        fw.close();
    }

    public enum MoveChoice {Real, NM, NMEC, Stage, StageEC};

    public static class SimParams extends ParameterBase {
        public int D = 3;
        public double k2 = 218.220183;
        public double gammaLangevin = Math.sqrt(k2);
        public long steps = 100000;
        public double density = 1.0;
        public double temperature = 1.0;
        public int numAtoms = 108;
        public double mass = 1.0;
        public double hbar = 0.1;
        public double rc = 2.5;
        public boolean isGraphic = false;
        public int nShifts = 0;
        public double timeStep = 0.01;
        public int nBeads = 1;
        public LJPIMD2.MoveChoice coordType = LJPIMD2.MoveChoice.StageEC;
        public double facTimestep = 0.01;
    }
}
