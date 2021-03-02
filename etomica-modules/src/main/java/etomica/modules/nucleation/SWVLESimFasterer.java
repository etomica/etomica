/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.nucleation;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMCFasterer;
import etomica.integrator.IntegratorManagerMC;
import etomica.integrator.mcmove.MCMoveAtomFasterer;
import etomica.integrator.mcmove.MCMoveMoleculeExchangeVLEFasterer;
import etomica.integrator.mcmove.MCMoveVolumeExchangeVLEFasterer;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.nbr.cell.PotentialMasterCellFasterer;
import etomica.potential.BondingInfo;
import etomica.potential.P2HardGeneric;
import etomica.potential.P2SquareWell;
import etomica.potential.PotentialMasterFasterer;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesGeneral;
import etomica.util.random.RandomMersenneTwister;

public class SWVLESimFasterer extends Simulation {

    public final Box boxLiquid, boxVapor;
    public final SpeciesGeneral species;
    public final IntegratorMCFasterer integratorLiquid, integratorVapor;
    public final IntegratorManagerMC integratorGEMC;

    protected final P2HardGeneric p2;
    protected double temperature;
    protected double density;

    public SWVLESimFasterer(int D) {
        super(Space.getInstance(D));
        setRandom(new RandomMersenneTwister(5));
        boolean doNBR = true;
        int initNumMolecules = 400;
        temperature = 1;
        density = 0.3;

        double initBoxSize = Math.pow(initNumMolecules / density, (1.0 / D));

        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species);

        boxLiquid = this.makeBox(new BoundaryRectangularPeriodic(space, initBoxSize));
        boxVapor = this.makeBox(new BoundaryRectangularPeriodic(space, initBoxSize));
        boxLiquid.setNMolecules(species, initNumMolecules);
        boxVapor.setNMolecules(species, initNumMolecules);
        Configuration config = new ConfigurationLattice(D == 2 ? new LatticeOrthorhombicHexagonal(space) : new LatticeCubicFcc(space), space);
        config.initializeCoordinates(boxLiquid);
        config.initializeCoordinates(boxVapor);

        final double range = 1.5;
        PotentialMasterFasterer potentialMasterV = new PotentialMasterFasterer(getSpeciesManager(), boxVapor, BondingInfo.noBonding());
        PotentialMasterFasterer potentialMasterL = new PotentialMasterFasterer(getSpeciesManager(), boxLiquid, BondingInfo.noBonding());
        p2 = P2SquareWell.makePotential(1.0, 1.5, 1.0);
        if (doNBR) {
            potentialMasterV = new PotentialMasterCellFasterer(getSpeciesManager(), boxVapor, 2, BondingInfo.noBonding());
            potentialMasterL = new PotentialMasterCellFasterer(getSpeciesManager(), boxLiquid, 2, BondingInfo.noBonding());
        } else {
            potentialMasterV = new PotentialMasterFasterer(getSpeciesManager(), boxVapor, BondingInfo.noBonding());
            potentialMasterL = new PotentialMasterFasterer(getSpeciesManager(), boxLiquid, BondingInfo.noBonding());
        }
        potentialMasterV.setPairPotential(species.getLeafType(), species.getLeafType(), p2);
        potentialMasterL.setPairPotential(species.getLeafType(), species.getLeafType(), p2);

        integratorLiquid = new IntegratorMCFasterer(potentialMasterL, random, temperature, boxLiquid);
        integratorLiquid.getMoveManager().setEquilibrating(true);
        MCMoveAtomFasterer atomMove = new MCMoveAtomFasterer(random, potentialMasterL, boxLiquid);
        integratorLiquid.getMoveManager().addMCMove(atomMove);

        integratorVapor = new IntegratorMCFasterer(potentialMasterV, random, temperature, boxVapor);
        integratorVapor.getMoveManager().setEquilibrating(true);
        atomMove = new MCMoveAtomFasterer(random, potentialMasterV, boxVapor);
        integratorVapor.getMoveManager().addMCMove(atomMove);

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
        final MCMoveVolumeExchangeVLEFasterer volumeExchange = new MCMoveVolumeExchangeVLEFasterer(random, space, integratorLiquid, integratorVapor);
        volumeExchange.setStepSize(0.05);
        MCMoveMoleculeExchangeVLEFasterer moleculeExchange = new MCMoveMoleculeExchangeVLEFasterer(
                random, space, integratorLiquid, integratorVapor);
        integratorGEMC.getMoveManager().addMCMove(volumeExchange);
        integratorGEMC.getMoveManager().addMCMove(moleculeExchange);
//        integratorGEMC.getMoveManager().setFrequency(volumeExchange, 0.01);

        getController().addActivity(new ActivityIntegrate(integratorGEMC));
    }

}
