/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.heVLE;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.Helium;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.IData;
import etomica.data.meter.MeterDensity;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorManagerMC;
import etomica.integrator.mcmove.*;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.NeighborCellManager;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.units.*;
import etomica.util.IListener;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

public class VLESim extends Simulation {

    public final Box boxLiquid, boxVapor;
    public final SpeciesGeneral species;
    public final IntegratorMC integratorLiquid, integratorVapor;
    public final IntegratorManagerMC integratorGEMC;
    protected final Potential2SoftSpherical p2;
    protected final P2SoftSphericalTruncated p2Truncated;
    protected double sigma;
    protected double temperature;
    protected double epsilon;
    protected double moment;
    protected double density;

    public VLESim(GEMCParams params) {
        super(Space3D.getInstance());
        boolean doNBR = false;
        int initNumMolecules = params.numAtoms;
        temperature = Kelvin.UNIT.toSim(params.temperatureK);
        CompoundUnit densityUnit = new CompoundUnit(new Unit[]{Mole.UNIT, Liter.UNIT}, new double[]{1, -1});
        density = densityUnit.toSim(params.density);

        double initBoxSize = Math.pow(initNumMolecules / density, (1.0 / 3.0));

        species = SpeciesGeneral.monatomic(space, AtomType.element(Helium.INSTANCE));
        addSpecies(species);

        System.out.println("box size: " + initBoxSize);

        boxLiquid = this.makeBox(new BoundaryRectangularPeriodic(space, initBoxSize));
        boxVapor = this.makeBox(new BoundaryRectangularPeriodic(space, initBoxSize));
        boxLiquid.setNMolecules(species, initNumMolecules);
        boxVapor.setNMolecules(species, initNumMolecules);
        Configuration config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        config.initializeCoordinates(boxLiquid);
        config.initializeCoordinates(boxVapor);

        final double range = 10.0;
        PotentialMaster potentialMaster = new PotentialMasterMonatomic(this);
        if (doNBR) {
            potentialMaster = new PotentialMasterCell(this, range, space);
            ((PotentialMasterCell) potentialMaster).setCellRange(2);
        }
        p2 = params.approx ? new P2HeSimplified(space) : new P2HePCKLJS(space);
        p2Truncated = new P2SoftSphericalTruncated(getSpace(), p2, range);
//        ((P2SoftSphericalTruncatedBox)potential).setTruncationFactor(0.35);
        potentialMaster.addPotential(p2Truncated, new AtomType[]{species.getLeafType(), species.getLeafType()});

        integratorLiquid = new IntegratorMC(potentialMaster, random, temperature, boxLiquid);
        integratorLiquid.getMoveManager().setEquilibrating(true);
        MCMoveAtom atomMove = new MCMoveAtom(random, potentialMaster, space, 0.5, 5.0, true);
        integratorLiquid.getMoveManager().addMCMove(atomMove);
//        ((MCMoveStepTracker)atomMove.getTracker()).setNoisyAdjustment(true);

        integratorVapor = new IntegratorMC(potentialMaster, random, temperature, boxVapor);
        integratorVapor.getMoveManager().setEquilibrating(true);
        atomMove = new MCMoveAtom(random, potentialMaster, space, 0.5, 10.0, true);
        integratorVapor.getMoveManager().addMCMove(atomMove);
//        ((MCMoveStepTracker)atomMove.getTracker()).setNoisyAdjustment(true);

        if (!doNBR) {
            BoxImposePbc pbc = new BoxImposePbc(boxLiquid, space);
            IntegratorListenerAction pbcListener = new IntegratorListenerAction(pbc);
            integratorLiquid.getEventManager().addListener(pbcListener);
            pbcListener.setInterval(100);
            pbc = new BoxImposePbc(boxVapor, space);
            pbcListener = new IntegratorListenerAction(pbc);
            integratorVapor.getEventManager().addListener(pbcListener);
            pbcListener.setInterval(100);
        }

        integratorGEMC = new IntegratorManagerMC(random);
        integratorGEMC.setTemperature(temperature);
        integratorGEMC.getMoveManager().setEquilibrating(true);
        integratorGEMC.setGlobalMoveInterval(5);
        integratorGEMC.addIntegrator(integratorLiquid);
        integratorGEMC.addIntegrator(integratorVapor);
        final MCMoveVolumeExchangeVLE volumeExchange = new MCMoveVolumeExchangeVLE(
                potentialMaster, random, space, integratorLiquid, integratorVapor);
        volumeExchange.setStepSize(0.05);
        MCMoveMoleculeExchangeVLE moleculeExchange = new MCMoveMoleculeExchangeVLE(
                potentialMaster, random, space, integratorLiquid, integratorVapor);
        integratorGEMC.getMoveManager().addMCMove(volumeExchange);
        integratorGEMC.getMoveManager().addMCMove(moleculeExchange);
        integratorGEMC.getMoveManager().setFrequency(volumeExchange, 0.01);

        integratorGEMC.getMoveEventManager().addListener(new IListener<MCMoveEvent>() {
            public void actionPerformed(MCMoveEvent event) {
                if (event instanceof MCMoveTrialCompletedEvent &&
                        ((MCMoveTrialCompletedEvent) event).isAccepted()) {
                    return;
                }
                if (event.getMCMove() == volumeExchange) {
                    if (boxLiquid.getBoundary().getBoxSize().getX(0) * 0.499 < range) {
                        p2Truncated.setTruncationRadius(0.499 * boxLiquid.getBoundary().getBoxSize().getX(0));
                    } else {
                        p2Truncated.setTruncationRadius(range);
                    }
                }
            }
        });

//        integratorGEMC.getEventManager().addListener(new IntegratorListener() {
//            boolean fixed = false;
//            public void integratorInitialized(IntegratorEvent e) {}
//            public void integratorStepStarted(IntegratorEvent e) {}
//            public void integratorStepFinished(IntegratorEvent e) {
//                if (!fixed && integratorGEMC.getStepCount() > 100000) {
//                    integratorGEMC.setGlobalMoveInterval(2);
//                }
//            }
//        });
//

        if (doNBR) {
            ((PotentialMasterCell) potentialMaster).getBoxCellManager(boxLiquid).assignCellAll();
            ((PotentialMasterCell) potentialMaster).getBoxCellManager(boxVapor).assignCellAll();
            integratorLiquid.getMoveEventManager().addListener(((NeighborCellManager) ((PotentialMasterCell) potentialMaster).getBoxCellManager(boxLiquid)).makeMCMoveListener());
            integratorVapor.getMoveEventManager().addListener((((NeighborCellManager) ((PotentialMasterCell) potentialMaster).getBoxCellManager(boxVapor)).makeMCMoveListener()));
        }
    }

    public static void main(String[] args) {
        GEMCParams params = new GEMCParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
        }
        final VLESim sim = new VLESim(params);

        long steps = params.numSteps;
        long t1 = System.currentTimeMillis();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorGEMC, steps / 10));
        sim.integratorGEMC.resetStepCount();

        int interval = params.numAtoms;
        int blockSize = (int) (params.numSteps / (interval * 100));
        MeterDensity liquidDensity = new MeterDensity(sim.boxLiquid);
        AccumulatorAverageFixed accLiquidDensity = new AccumulatorAverageFixed(blockSize);
        DataPumpListener pumpLiquidDensity = new DataPumpListener(liquidDensity, accLiquidDensity, interval);
        sim.integratorLiquid.getEventManager().addListener(pumpLiquidDensity);

        MeterDensity vaporDensity = new MeterDensity(sim.boxVapor);
        AccumulatorAverageFixed accVaporDensity = new AccumulatorAverageFixed(blockSize);
        DataPumpListener pumpVaporDensity = new DataPumpListener(vaporDensity, accVaporDensity, interval);
        sim.integratorLiquid.getEventManager().addListener(pumpVaporDensity);

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorGEMC, steps));
        long t2 = System.currentTimeMillis();

        IData liquidDensityData = accLiquidDensity.getData();
        double liquidAvg = liquidDensityData.getValue(accLiquidDensity.AVERAGE.index);
        double liquidErr = liquidDensityData.getValue(accLiquidDensity.ERROR.index);
        double liquidCor = liquidDensityData.getValue(accLiquidDensity.BLOCK_CORRELATION.index);

        System.out.println("liquid density: " + liquidAvg + " err: " + liquidErr + " cor: " + liquidCor);

        IData vaporDensityData = accVaporDensity.getData();
        double vaporAvg = vaporDensityData.getValue(accVaporDensity.AVERAGE.index);
        double vaporErr = vaporDensityData.getValue(accVaporDensity.ERROR.index);
        double vaporCor = vaporDensityData.getValue(accVaporDensity.BLOCK_CORRELATION.index);
        System.out.println("vapor density: " + vaporAvg + " err: " + vaporErr + " cor: " + vaporCor);

        System.out.println("time: " + (t2 - t1) * 0.001);
    }

    public static class GEMCParams extends ParameterBase {
        public double temperatureK = 10;
        public int numAtoms = 200;
        public double density = 25; // kmol/m^3 = mol/L
        public boolean approx = true;
        public long numSteps = 1000000;
    }
}
