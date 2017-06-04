/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.vle;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.IAtomType;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorManagerMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveEvent;
import etomica.integrator.mcmove.MCMoveRotate;
import etomica.integrator.mcmove.MCMoveTrialCompletedEvent;
import etomica.lattice.LatticeCubicFcc;
import etomica.listener.IntegratorListenerAction;
import etomica.nbr.cell.NeighborCellManager;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2LJQ;
import etomica.potential.P2SoftTruncated;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresRotating;
import etomica.units.Debye;
import etomica.units.Kelvin;
import etomica.util.IEvent;
import etomica.util.IListener;

public class VLESim extends Simulation {

    public final Box boxLiquid, boxVapor;
    public final SpeciesSpheresRotating species;
    public final IntegratorMC integratorLiquid, integratorVapor;
    public final IntegratorManagerMC integratorGEMC;
    public final ActivityIntegrate activityIntegrate;
    protected double sigma;
    protected double temperature;
    protected double epsilon;
    protected double moment;
    protected double density;
    protected final P2LJQ p2LJQ;
    protected final P2SoftTruncated p2Truncated;
    
    public static void main(String[] args) {
        new VLESim();
    }
    
    public VLESim() {
        super(Space3D.getInstance());
        boolean doNBR = false;
        int initNumMolecules = 200;
        sigma = 3;
        temperature = Kelvin.UNIT.toSim(250);
        epsilon = Kelvin.UNIT.toSim(150);
        moment = Debye.UNIT.toSim(5);
        moment *= moment;
        density = 0.004;

        double initBoxSize = Math.pow(initNumMolecules/density, (1.0/3.0));
        
        species = new SpeciesSpheresRotating(this, space);
        addSpecies(species);

        boxLiquid = new Box(new BoundaryRectangularPeriodic(space, initBoxSize), space);
        addBox(boxLiquid);
        boxVapor = new Box(new BoundaryRectangularPeriodic(space, initBoxSize), space);
        addBox(boxVapor);
        boxLiquid.setNMolecules(species, initNumMolecules);
        boxVapor.setNMolecules(species, initNumMolecules);
        Configuration config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        config.initializeCoordinates(boxLiquid);
        config.initializeCoordinates(boxVapor);

        final double range = 15.0;
        PotentialMaster potentialMaster = new PotentialMasterMonatomic(this);
        if (doNBR) {
            potentialMaster = new PotentialMasterCell(this, range, space);
            ((PotentialMasterCell)potentialMaster).setCellRange(2);
        }
        p2LJQ = new P2LJQ(space, sigma, epsilon, moment);
        p2LJQ.setTemperature(temperature);
        p2Truncated = new P2SoftTruncated(p2LJQ, range, space);
//        ((P2SoftSphericalTruncatedBox)potential).setTruncationFactor(0.35);
        potentialMaster.addPotential(p2Truncated, new IAtomType[]{species.getLeafType(), species.getLeafType()});
        
        integratorLiquid = new IntegratorMC(potentialMaster, random, temperature);
        integratorLiquid.getMoveManager().setEquilibrating(true);
        integratorLiquid.setBox(boxLiquid);
        MCMoveAtom atomMove = new MCMoveAtom(potentialMaster, random, space, 0.5, 5.0, true);
        integratorLiquid.getMoveManager().addMCMove(atomMove);
        MCMoveRotate rotateMove = new MCMoveRotate(potentialMaster, random, space);
        integratorLiquid.getMoveManager().addMCMove(rotateMove);
//        ((MCMoveStepTracker)atomMove.getTracker()).setNoisyAdjustment(true);
        
        integratorVapor = new IntegratorMC(potentialMaster, random, temperature);
        integratorVapor.setBox(boxVapor);
        integratorVapor.getMoveManager().setEquilibrating(true);
        atomMove = new MCMoveAtom(potentialMaster, random, space, 0.5, 5.0, true);
        integratorVapor.getMoveManager().addMCMove(atomMove);
        rotateMove = new MCMoveRotate(potentialMaster, random, space);
        integratorVapor.getMoveManager().addMCMove(rotateMove);
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
        integratorGEMC.getMoveManager().setEquilibrating(true);
        integratorGEMC.setGlobalMoveInterval(2);
        integratorGEMC.addIntegrator(integratorLiquid);
        integratorGEMC.addIntegrator(integratorVapor);
        final MCMoveVolumeExchangeVLE volumeExchange = new MCMoveVolumeExchangeVLE(
                potentialMaster, random, space, integratorLiquid,integratorVapor);
        volumeExchange.setStepSize(0.05);
        MCMoveMoleculeExchangeVLE moleculeExchange = new MCMoveMoleculeExchangeVLE(
                potentialMaster, random, space, integratorLiquid,integratorVapor);
        integratorGEMC.getMoveManager().addMCMove(volumeExchange);
        integratorGEMC.getMoveManager().addMCMove(moleculeExchange);
        integratorGEMC.getMoveManager().setFrequency(volumeExchange, 0.01);

        integratorGEMC.getMoveEventManager().addListener(new IListener() {
            public void actionPerformed(IEvent event) {
                if (event instanceof MCMoveTrialCompletedEvent &&
                    ((MCMoveTrialCompletedEvent)event).isAccepted()) {
                    return;
                }
                if (((MCMoveEvent)event).getMCMove() == volumeExchange) {
                    if (boxLiquid.getBoundary().getBoxSize().getX(0)*0.499 < range) {
                        p2Truncated.setTruncationRadius(0.499*boxLiquid.getBoundary().getBoxSize().getX(0));
                    }
                    else {
                        p2Truncated.setTruncationRadius(range);
                    }
                }
            }
        });

        activityIntegrate = new ActivityIntegrate(integratorGEMC);
        getController().addAction(activityIntegrate);

        if (doNBR) {
            ((PotentialMasterCell)potentialMaster).getCellAgentManager().getAgent(boxLiquid).assignCellAll();
            ((PotentialMasterCell)potentialMaster).getCellAgentManager().getAgent(boxVapor).assignCellAll();
            integratorLiquid.getMoveEventManager().addListener(((NeighborCellManager)((PotentialMasterCell)potentialMaster).getCellAgentManager().getAgent(boxLiquid)).makeMCMoveListener());
            integratorVapor.getMoveEventManager().addListener(((NeighborCellManager)((PotentialMasterCell)potentialMaster).getCellAgentManager().getAgent(boxVapor)).makeMCMoveListener());
        }
    }
    
    public void setSigma(double newSigma) {
        sigma = newSigma;
        p2LJQ.setSigma(sigma);
        p2Truncated.setTruncationRadius(4.0*sigma);
        integratorLiquid.reset();
        integratorVapor.reset();
    }

    public double getSigma() {
        return sigma;
    }
    
    public void setEpsilon(double newEpsilon) {
        p2LJQ.setEpsilon(newEpsilon);
        integratorLiquid.reset();
        integratorVapor.reset();
    }
    
    public double getEpsilon() {
        return p2LJQ.getEpsilon();
    }
    
    public void setMoment(double newQ) {
        p2LJQ.setQuadrupolarMomentSquare(newQ*newQ);
        integratorLiquid.reset();
        integratorVapor.reset();
    }
    
    public double getMoment() {
        return Math.sqrt(p2LJQ.getQuadrupolarMomentSquare());
    }
}
