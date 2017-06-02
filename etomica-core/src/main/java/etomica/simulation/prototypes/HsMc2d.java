/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation.prototypes;
import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.api.IAtomType;
import etomica.box.Box;
import etomica.potential.PotentialMaster;
import etomica.config.ConfigurationLattice;
import etomica.data.DataSourceCountSteps;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.listener.IntegratorListenerAction;
import etomica.potential.P2HardSphere;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space2d.Space2D;
import etomica.species.SpeciesSpheresMono;

/**
 * Simple hard-sphere Monte Carlo simulation in 2D.
 *
 * @author David Kofke
 */
 
public class HsMc2d extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public IntegratorMC integrator;
    public MCMoveAtom mcMoveAtom;
    public SpeciesSpheresMono species, species2;
    public Box box;
    public P2HardSphere potential;
    public Controller controller;
    public DataSourceCountSteps meterCycles;

    public HsMc2d() {
        super(Space2D.getInstance());
        PotentialMaster potentialMaster = new PotentialMasterMonatomic(this);
	    integrator = new IntegratorMC(this, potentialMaster);
	    mcMoveAtom = new MCMoveAtom(random, potentialMaster, space);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        species = new SpeciesSpheresMono(this, space);
        species2 = new SpeciesSpheresMono(this, space);
        addSpecies(species);
        addSpecies(species2);
        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, 20);
        box.setNMolecules(species2, 20);
        new ConfigurationLattice(new LatticeOrthorhombicHexagonal(space), space).initializeCoordinates(box);
	    potential = new P2HardSphere(space);
	    IAtomType type1 = species.getLeafType();
        IAtomType type2 = species2.getLeafType();
        potentialMaster.addPotential(potential, new IAtomType[] {type1, type1});
        potentialMaster.addPotential(potential, new IAtomType[] {type1, type2});
        potentialMaster.addPotential(potential, new IAtomType[] {type2, type2});
	    meterCycles = new DataSourceCountSteps(integrator);

        integrator.setBox(box);
        integrator.getMoveManager().addMCMove(mcMoveAtom);
        integrator.getEventManager().addListener(new IntegratorListenerAction(new BoxImposePbc(box, space)));

//	    LatticeRenderer.ColorSchemeCell colorSchemeCell = new LatticeRenderer.ColorSchemeCell();
//	    display.setColorScheme(colorSchemeCell);
	    
//		elementCoordinator.go();
//	    etomica.lattice.BravaisLattice lattice = ((IteratorFactoryCell)this.getIteratorFactory()).getLattice(box);
//        colorSchemeCell.setLattice(lattice);
    }
    
}
