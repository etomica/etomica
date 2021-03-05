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
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorManagerMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveMoleculeExchangeVLE;
import etomica.integrator.mcmove.MCMoveVolumeExchangeVLE;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.nbr.cell.NeighborCellManager;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2SquareWell;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesGeneral;
import etomica.util.random.RandomMersenneTwister;

public class SWVLESim extends Simulation {

    public final Box boxLiquid, boxVapor;
    public final SpeciesGeneral species;
    public final IntegratorMC integratorLiquid, integratorVapor;
    public final IntegratorManagerMC integratorGEMC;

    protected final P2SquareWell p2;
    protected double temperature;
    protected double density;

    public SWVLESim(int D) {
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
        PotentialMaster potentialMaster = new PotentialMasterMonatomic(getSpeciesManager());
        if (doNBR) {
            potentialMaster = new PotentialMasterCell(this, range);
            ((PotentialMasterCell) potentialMaster).setCellRange(2);
        }
        p2 = new P2SquareWell(space, 1.0, 1.5, 1.0, false);
        potentialMaster.addPotential(p2, new AtomType[]{species.getLeafType(), species.getLeafType()});

        integratorLiquid = new IntegratorMC(potentialMaster, random, temperature, boxLiquid);
        integratorLiquid.getMoveManager().setEquilibrating(true);
        MCMoveAtom atomMove = new MCMoveAtom(random, potentialMaster, space, 0.5, 5.0, true);
        integratorLiquid.getMoveManager().addMCMove(atomMove);

        integratorVapor = new IntegratorMC(potentialMaster, random, temperature, boxVapor);
        integratorVapor.getMoveManager().setEquilibrating(true);
        atomMove = new MCMoveAtom(random, potentialMaster, space, 0.5, 5.0, true);
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
        final MCMoveVolumeExchangeVLE volumeExchange = new MCMoveVolumeExchangeVLE(
                potentialMaster, random, space, integratorLiquid, integratorVapor);
        volumeExchange.setStepSize(0.05);
        MCMoveMoleculeExchangeVLE moleculeExchange = new MCMoveMoleculeExchangeVLE(
                potentialMaster, random, space, integratorLiquid, integratorVapor);
        integratorGEMC.getMoveManager().addMCMove(volumeExchange);
        integratorGEMC.getMoveManager().addMCMove(moleculeExchange);
//        integratorGEMC.getMoveManager().setFrequency(volumeExchange, 0.01);

        getController().addActivity(new ActivityIntegrate(integratorGEMC));

        if (doNBR) {
            ((PotentialMasterCell) potentialMaster).getBoxCellManager(boxLiquid).assignCellAll();
            ((PotentialMasterCell) potentialMaster).getBoxCellManager(boxVapor).assignCellAll();
            integratorLiquid.getMoveEventManager().addListener(((NeighborCellManager) ((PotentialMasterCell) potentialMaster).getBoxCellManager(boxLiquid)).makeMCMoveListener());
            integratorVapor.getMoveEventManager().addListener((((NeighborCellManager) ((PotentialMasterCell) potentialMaster).getBoxCellManager(boxVapor)).makeMCMoveListener()));
        }
    }

}
