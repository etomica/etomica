/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.vle;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMCFasterer;
import etomica.integrator.IntegratorManagerMC;
import etomica.integrator.mcmove.*;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.NeighborCellManagerFasterer;
import etomica.potential.BondingInfo;
import etomica.potential.P2LJQ;
import etomica.potential.P2SoftTruncated;
import etomica.potential.compute.NeighborManager;
import etomica.potential.compute.NeighborManagerSimple;
import etomica.potential.compute.PotentialComputePairGeneral;
import etomica.simulation.Simulation;
import etomica.space.BoundaryEvent;
import etomica.space.BoundaryEventListener;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.species.SpeciesSpheresRotating;
import etomica.units.Debye;
import etomica.units.Kelvin;
import etomica.util.random.RandomMersenneTwister;

public class VLESimFasterer extends Simulation {

    public final Box boxLiquid, boxVapor;
    public final SpeciesGeneral species;
    public final IntegratorMCFasterer integratorLiquid, integratorVapor;
    public final IntegratorManagerMC integratorGEMC;
    protected final P2LJQ p2LJQ;
    protected final P2SoftTruncated p2Truncated;
    protected double sigma;
    protected double temperature;
    protected double epsilon;
    protected double moment;
    protected double density;

    public VLESimFasterer() {
        super(Space3D.getInstance());
        setRandom(new RandomMersenneTwister(2));
        boolean doNBR = false;
        int initNumMolecules = 200;
        sigma = 3;
        temperature = Kelvin.UNIT.toSim(250);
        epsilon = Kelvin.UNIT.toSim(150);
        moment = Debye.UNIT.toSim(5);
        moment *= moment;
        density = 0.004;

        double initBoxSize = Math.pow(initNumMolecules / density, (1.0 / 3.0));

        species = SpeciesSpheresRotating.create(space, new ElementSimple(this), false, true);
        addSpecies(species);

        boxLiquid = this.makeBox(new BoundaryRectangularPeriodic(space, initBoxSize));
        boxVapor = this.makeBox(new BoundaryRectangularPeriodic(space, initBoxSize));
        boxLiquid.setNMolecules(species, initNumMolecules);
        boxVapor.setNMolecules(species, initNumMolecules);
        Configuration config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        config.initializeCoordinates(boxLiquid);
        config.initializeCoordinates(boxVapor);

        final double range = 4.0 * sigma;
        NeighborManager neighborManagerL = doNBR ?
                new NeighborCellManagerFasterer(getSpeciesManager(), boxLiquid, 2, BondingInfo.noBonding())
                : new NeighborManagerSimple(boxLiquid);
        NeighborManager neighborManagerV = doNBR ?
                new NeighborCellManagerFasterer(getSpeciesManager(), boxVapor, 2, BondingInfo.noBonding())
                : new NeighborManagerSimple(boxVapor);
        PotentialComputePairGeneral potentialMasterL = new PotentialComputePairGeneral(getSpeciesManager(), boxLiquid, neighborManagerL);
        PotentialComputePairGeneral potentialMasterV = new PotentialComputePairGeneral(getSpeciesManager(), boxVapor, neighborManagerV);
        p2LJQ = new P2LJQ(space, sigma, epsilon, moment);
        p2LJQ.setTemperature(temperature);
        p2Truncated = new P2SoftTruncated(p2LJQ, range, space);
//        ((P2SoftSphericalTruncatedBox)potential).setTruncationFactor(0.35);
        potentialMasterL.setPairPotential(species.getLeafType(), species.getLeafType(), p2Truncated);
        potentialMasterV.setPairPotential(species.getLeafType(), species.getLeafType(), p2Truncated);

        integratorLiquid = new IntegratorMCFasterer(potentialMasterL, random, temperature, boxLiquid);
        integratorLiquid.getMoveManager().setEquilibrating(true);
        MCMoveAtomFasterer atomMoveL = new MCMoveAtomFasterer(random, potentialMasterL, boxLiquid);//space, 0.5, 5.0, true);
        atomMoveL.setStepSize(0.5);
        atomMoveL.setStepSizeMax(5.0);
        integratorLiquid.getMoveManager().addMCMove(atomMoveL);
        MCMoveAtomRotateFasterer rotateMoveL = new MCMoveAtomRotateFasterer(random, potentialMasterL, boxLiquid);
        integratorLiquid.getMoveManager().addMCMove(rotateMoveL);
        ((MCMoveStepTracker) atomMoveL.getTracker()).setNoisyAdjustment(true);

        integratorVapor = new IntegratorMCFasterer(potentialMasterV, random, temperature, boxVapor);
        integratorVapor.getMoveManager().setEquilibrating(true);
        MCMoveAtomFasterer atomMoveV = new MCMoveAtomFasterer(random, potentialMasterV, boxLiquid);
        atomMoveV.setStepSize(1.0);
        atomMoveV.setStepSizeMax(10.0);
        integratorVapor.getMoveManager().addMCMove(atomMoveV);
        MCMoveAtomRotateFasterer rotateMoveV = new MCMoveAtomRotateFasterer(random, potentialMasterV, boxVapor);
        integratorVapor.getMoveManager().addMCMove(rotateMoveV);
        ((MCMoveStepTracker) atomMoveV.getTracker()).setNoisyAdjustment(true);

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
        integratorGEMC.setGlobalMoveInterval(2);
        integratorGEMC.addIntegrator(integratorLiquid);
        integratorGEMC.addIntegrator(integratorVapor);
        final MCMoveVolumeExchangeVLEFasterer volumeExchange = new MCMoveVolumeExchangeVLEFasterer(
                random, space, integratorLiquid, integratorVapor);
        volumeExchange.setStepSize(0.05);
        MCMoveMoleculeExchangeFasterer moleculeExchange = new MCMoveMoleculeExchangeFasterer(
                random, space, integratorLiquid, integratorVapor);
        integratorGEMC.getMoveManager().addMCMove(volumeExchange);
        integratorGEMC.getMoveManager().addMCMove(moleculeExchange);
        integratorGEMC.getMoveManager().setFrequency(volumeExchange, 0.01);


        boxLiquid.getBoundary().getEventManager().addListener(new BoundaryEventListener() {
            @Override
            public void boundaryInflate(BoundaryEvent e) {
                double rc = Math.min(4.0 * sigma, boxLiquid.getBoundary().getBoxSize().getX(0) * 0.499);
                p2Truncated.setTruncationRadius(rc);
            }
        });
        this.getController().addActivity(new ActivityIntegrate(integratorGEMC));
    }

    public static void main(String[] args) {
        new VLESimFasterer();
    }

    public double getSigma() {
        return sigma;
    }

    public void setSigma(double newSigma) {
        sigma = newSigma;
        p2LJQ.setSigma(sigma);
        double rc = Math.min(4.0 * sigma, boxLiquid.getBoundary().getBoxSize().getX(0) * 0.499);
        p2Truncated.setTruncationRadius(rc);
        integratorLiquid.reset();
        integratorVapor.reset();
    }

    public double getEpsilon() {
        return p2LJQ.getEpsilon();
    }

    public void setEpsilon(double newEpsilon) {
        p2LJQ.setEpsilon(newEpsilon);
        integratorLiquid.reset();
        integratorVapor.reset();
    }

    public double getMoment() {
        return Math.sqrt(p2LJQ.getQuadrupolarMomentSquare());
    }

    public void setMoment(double newQ) {
        p2LJQ.setQuadrupolarMomentSquare(newQ * newQ);
        integratorLiquid.reset();
        integratorVapor.reset();
    }
}
