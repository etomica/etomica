/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.*;
import etomica.box.Box;
import etomica.config.ConformationLinear;
import etomica.data.AccumulatorAverageCovariance;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceScalar;
import etomica.data.IData;
import etomica.data.types.DataGroup;
import etomica.graphics.*;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorLangevin;
import etomica.integrator.IntegratorMD;
import etomica.potential.P1Anharmonic234;
import etomica.potential.P1AnharmonicTIA;
import etomica.potential.P2Harmonic;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputeField;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space1d.Space1D;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;
import etomica.units.dimensions.Length;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;
import java.util.ArrayList;
import java.util.List;

public class SimQuantumAOPIMD2 extends Simulation {

    public PotentialComputeField pcP1, pcP2;
    public final IntegratorMD integrator;
    public final Box box;
    public P1Anharmonic234 p1ah;
    public double betaN;
    public int nBeads;
    public MCMoveHOReal2 moveStageSimple, moveStageEC;
    public int dim;

    public SimQuantumAOPIMD2(Space space, MoveChoice coordType, double mass, double timeStep, double gammaLangevin, int nBeads, double temperature, double k3, double k4, double omega, boolean isTIA, double hbar, boolean isExactA, boolean isCayleyA) {
        super(space);

        if (isExactA && isCayleyA) {
            System.out.println(" Can not have both isExactA and isCayleyA to be true!");
            System.exit(0);
        }
        SpeciesGeneral species = new SpeciesBuilder(space)
                .setDynamic(true)
                .addCount(AtomType.simple("A", mass / nBeads), nBeads)
                .withConformation(new ConformationLinear(space, 0))
                .build();
        addSpecies(species);
        SpeciesGeneral speciesLattice = null;
        if (space.D() > 1) {
            speciesLattice = new SpeciesBuilder(space)
                    .setDynamic(true)
                    .addAtom(AtomType.simple("L", Double.POSITIVE_INFINITY), space.makeVector())
                    .build();
            addSpecies(speciesLattice);
        }
        this.nBeads = nBeads;
        box = this.makeBox(new BoundaryRectangularNonperiodic(space));
        box.setNMolecules(species, 1);
        if (speciesLattice!=null) box.setNMolecules(speciesLattice, 1);
        //pm2 that uses the full PI potential, for data collection
        //spring P2 part (x_i-x_{i+1})^2
        this.dim = space.D();
        double beta = 1.0/temperature;
        betaN = beta/nBeads;
        double omega2 = omega*omega;

        pcP1 = new PotentialComputeField(getSpeciesManager(), box);
        pcP2 = new PotentialComputeField(getSpeciesManager(), box);



        p1ah = new P1Anharmonic234(space, omega2/nBeads, k3/nBeads, k4/nBeads);
        pcP1.setFieldPotential(species.getLeafType(), p1ah);

        PotentialComputeAggregate.localStorageDefault = true;

        if (coordType == MoveChoice.Real) {
            integrator = new IntegratorLangevin(pcP1, random, timeStep, temperature, box, gammaLangevin);
        } else if (coordType == MoveChoice.NM) {
            MCMoveHO move = new MCMoveHO(space, pcP1, random, temperature, 0, box, hbar);
            if (isExactA) {
                integrator = new IntegratorLangevinPINMExactA(pcP1, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
            } else if (isCayleyA) {
                integrator = new IntegratorLangevinPINMCayleyA(pcP1, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
            } else {
                integrator = new IntegratorLangevinPINM(pcP1, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
            }
        } else if (coordType == MoveChoice.NMEC) {
            MCMoveHO move = new MCMoveHO(space, pcP1, random, temperature, omega2, box, hbar);
            if (isExactA) {
                integrator = new IntegratorLangevinPINMExactA(pcP2, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
            } else if (isCayleyA) {
                integrator = new IntegratorLangevinPINMCayleyA(pcP2, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
            } else {
                integrator = new IntegratorLangevinPINM(pcP1, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
            }
        } else if (coordType == MoveChoice.Stage) {
            MCMoveHOReal2 move = new MCMoveHOReal2(space, pcP1, random, temperature, 0, box, hbar);
            if (isExactA) {
                integrator = new IntegratorLangevinPIExactA(pcP1, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
            } else if (isCayleyA) {
                integrator = new IntegratorLangevinPICayleyA(pcP1, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
            } else {
                integrator = new IntegratorLangevinPI(pcP1, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
            }
        } else { //StageEC -- default
            MCMoveHOReal2 move = new MCMoveHOReal2(space, pcP1, random, temperature, omega2, box, hbar);
            if (isExactA) {
                integrator = new IntegratorLangevinPIExactA(pcP2, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
            } else if (isCayleyA) {
                integrator = new IntegratorLangevinPICayleyA(pcP2, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
            } else {
                integrator = new IntegratorLangevinPI(pcP1, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
            }
        }

        moveStageSimple = new MCMoveHOReal2(space, pcP1, random, temperature, 0, box, hbar);
        moveStageEC = new MCMoveHOReal2(space, pcP1, random, temperature, omega2, box, hbar);
        integrator.setThermostatNoDrift(false);
        integrator.setIsothermal(true);
    }

    public Integrator getIntegrator() {
        return integrator;
    }


    public static void main(String[] args) {
        final long startTime = System.currentTimeMillis();
        OctaneParams params = new OctaneParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            // custom parameters
            params.hbar = 1.0;
            params.steps = 1000000;
            params.temperature = 1;
            params.omega = 1;
            params.k3 = 0;
            params.k4 = 0;
            params.gammaLangevin = params.omega;
            params.coordType = MoveChoice.StageEC;
            params.timeStep = 0.1;
            params.nBeads = 1;
            params.isGraphic = !true;
        }

        int nShifts = params.nShifts;
        double mass = params.mass;
        double temperature = params.temperature;
        double hbar = params.hbar;
        double omega = params.omega;
        double k3 = params.k3;
        double k4 = params.k4;
        double gammaLangevin = params.gammaLangevin;
        boolean isGraphic = params.isGraphic;
        boolean isExactA = params.isExactA;
        boolean isCayley = params.isCayleyA;
        long steps = params.steps;
        long stepsEq = steps/10;
        boolean isTIA = params.isTIA;
        MoveChoice coordType = params.coordType;
        double omega2 = omega*omega;


        final SimQuantumAOPIMD2 sim = new SimQuantumAOPIMD2(Space1D.getInstance(), coordType, mass, params.timeStep, gammaLangevin, params.nBeads, temperature, k3, k4, omega, isTIA, hbar, isExactA, isCayley);
        sim.integrator.reset();

        System.out.println(" PIMD-" + coordType);
        System.out.println(" mass: " + mass);
        System.out.println(" T: " + temperature);
        System.out.println(" hbar: " + hbar);
        System.out.println(" w: " + Math.sqrt(omega2));
        System.out.println(" x = beta*hbar*w = " + hbar*omega/temperature);
        System.out.println(" nShifts: "+ nShifts);
        System.out.println(" steps: " +  steps + " stepsEq: " + stepsEq);
        System.out.println(" timestep: " + params.timeStep);
        System.out.println(" k3: " + k3);
        System.out.println(" k4: " + k4);
        System.out.println(" isTIA: " + isTIA);
        System.out.println(" gammaLangevin: " + gammaLangevin);
        System.out.println(" isExactA: " + isExactA);
        System.out.println(" isCayley: " + isCayley);

        MeterPIVir2 meterVir = meterVir = new MeterPIVir2(sim.pcP1, temperature, sim.box);


        if (isGraphic) {
            sim.getController().addActivity(new ActivityIntegrate(sim.integrator));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.GRAPHIC_ONLY);
            int intervalG = 100;
            simGraphic.setPaintInterval(sim.box, intervalG);
            int nBeads = 1;
            int finalNBeads = nBeads;
            ColorScheme colorScheme = new ColorScheme() {
                protected Color[] allColors;

                public Color getAtomColor(IAtom a) {
                    if (a.getType().getMass() == Double.POSITIVE_INFINITY) return sim.space.D() == 2 ? Color.BLACK : Color.WHITE;
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

            DisplayBox displayBox = simGraphic.getDisplayBox(sim.box);
            displayBox.setColorScheme(colorScheme);
            ((DiameterHashByType) displayBox.getDiameterHash()).setDiameter(sim.getSpecies(0).getAtomType(0), 0.6);
//            ((DiameterHashByType) displayBox.getDiameterHash()).setDiameter(sim.getSpecies(1).getAtomType(0), 0.1);

            if (sim.space.D() == 3) {
                AtomPair pair = new AtomPair();
                for (int j = 0; j < 1; j++) {
                    IAtomList beads = sim.box.getMoleculeList().get(j).getChildList();
                    for (int i = 0; i < nBeads; i++) {
                        pair.atom0 = beads.get(i);
                        int next = i + 1;
                        if (next == nBeads) next = 0;
                        pair.atom1 = beads.get(next);
                        ((DisplayBoxCanvasG3DSys) displayBox.canvas).makeBond(pair, null);
                    }
                }
                IAtomList beads = sim.box.getLeafList();
                for (int i = 0; i < nBeads; i++) {
                    pair.atom0 = beads.get(i);
                    pair.atom1 = beads.get(nBeads);
                    ((DisplayBoxCanvasG3DSys) displayBox.canvas).makeBond(pair, Color.BLUE);
                }
            }
            else if (sim.space.D() == 2) {
                for (int j = 0; j < 1; j++) {
                    IAtomList beads = sim.box.getMoleculeList().get(j).getChildList();
                    for (int i = 0; i < nBeads; i++) {
                        int next = i + 1;
                        if (next == nBeads) next = 0;
                        displayBox.addDrawable(new SimQuantumAOPIMD.MyBond(beads.get(i).getPosition(), beads.get(next).getPosition(), Color.RED, sim.box.getBoundary()));
                    }
                }
                IAtomList beads = sim.box.getLeafList();
                for (int i = 0; i < nBeads; i++) {
                    displayBox.addDrawable(new SimQuantumAOPIMD.MyBond(beads.get(i).getPosition(), beads.get(nBeads).getPosition(), Color.BLUE, sim.box.getBoundary()));
                }
            }

            simGraphic.makeAndDisplayFrame("PIMD - "+coordType);

            DataSourceScalar meterCOM = new DataSourceScalar("COM", Length.DIMENSION) {
                @Override
                public double getDataAsScalar() {
                    Vector COM = sim.box.getSpace().makeVector();
                    for (IAtom atom : sim.box.getLeafList()) {
                        COM.PE(atom.getPosition());
                    }
                    COM.TE(1.0/sim.box.getLeafList().size());
                    return COM.getX(0);
                }
            };

            return;
        }
        System.out.flush();
        // equilibration
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, stepsEq));
        System.out.println("\n equilibration finished");
        int interval = 1;
        int blocks = 100;
        long blockSize = steps / (interval * blocks);
        if (blockSize == 0) blockSize = 1;
        System.out.println(" numBlocks: " + blocks + " blocksize: " + blockSize + " interval: " + interval);


        AccumulatorAverageCovariance accumulatorVir = new AccumulatorAverageCovariance(blockSize);
        if (meterVir != null) {
            DataPumpListener accumulatorPumpVir = new DataPumpListener(meterVir, accumulatorVir, interval);
            sim.integrator.getEventManager().addListener(accumulatorPumpVir);
        }

        //Short Run
        long stepsShort = steps/10;
        sim.integrator.resetStepCount();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, stepsShort));
        System.out.println(" End of short run of " + stepsShort + " steps");

        double beta = 1/temperature;
        DataGroup dataVir = (DataGroup) accumulatorVir.getData();
        IData dataAvgVir = dataVir.getData(accumulatorVir.AVERAGE.index);
        IData dataErrVir = dataVir.getData(accumulatorVir.ERROR.index);
        IData dataCovVir = dataVir.getData(accumulatorVir.COVARIANCE.index);

        double xHx = dataAvgVir.getValue(1);
        double sdot = -1.0/beta/(1.0+beta*xHx);
        System.out.println(" sdot: " + sdot);
        meterVir.setSdot(sdot);
        accumulatorVir.reset();

        // Run
        sim.integrator.resetStepCount();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps));

        // Vir

        double FX2 = dataAvgVir.getValue(2);
        double varFX = dataCovVir.getValue(0);
        double uconv = dataAvgVir.getValue(3);
        double uhma = dataAvgVir.getValue(4);
        double uhma2 = dataAvgVir.getValue(5);
        double uhma3 = dataAvgVir.getValue(6);
        double errUconv = dataErrVir.getValue(3);
        double errUhma = dataErrVir.getValue(4);
        double errUhma2 = dataErrVir.getValue(5);
        double errUhma3 = dataErrVir.getValue(6);

        System.out.println(" fx2: " + (-1.0/beta/(beta*beta*FX2-1)));
        System.out.println(" var: " + (-1.0/beta/(beta*beta*varFX)));
        System.out.println();

        System.out.println(" uConv: " + uconv + " " + errUconv);
        System.out.println(" uHMA: " + uhma + " " + errUhma);
        System.out.println(" uHMA2: " + uhma2 + " " + errUhma2);
        System.out.println(" uHMA3: " + uhma3 + " " + errUhma3);

        long endTime = System.currentTimeMillis();
        System.out.println("\n time: (min) " + (endTime - startTime)/60.0/1000.0);
    }

    public enum MoveChoice {Real, NM, NMEC, Stage, StageEC};

    public static class OctaneParams extends ParameterBase {
        public double temperature = 1;
        public double hbar = 1;
        public double omega = 1;
        public double gammaLangevin = omega;
        public double k3 = 0;
        public double k4 = 0;
        public long steps = 10_000_000;
        public boolean isGraphic = false;
        public boolean isTIA = false;
        public double mass = 1.0;
        public MoveChoice coordType = MoveChoice.Real;
        public double timeStep = 0.1;
        public int nBeads = 1;
        public int nShifts = 0;
        public boolean isExactA = false;
        public boolean isCayleyA = false;
    }
}