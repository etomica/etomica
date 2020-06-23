/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.insertion;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMD.ThermostatType;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.random.RandomMersenneTwister;

public class Insertion extends Simulation {
    
    public SpeciesSpheresMono species, speciesGhost;
    public Box box;
    public IntegratorHard integrator;
    public P2HardWrapper potentialWrapper;
    public P2DoubleWell potentialGhost;
    public ActivityIntegrate activityIntegrate;
    
    public Insertion(Space _space) {
        super(_space);
        setRandom(new RandomMersenneTwister(2));

        // species
        species = new SpeciesSpheresMono(this, space);//index 1
        species.setIsDynamic(true);
        addSpecies(species);
        speciesGhost = new SpeciesSpheresMono(this, space);
        speciesGhost.setIsDynamic(true);
        ((ElementSimple) speciesGhost.getLeafType().getElement()).setMass(Double.POSITIVE_INFINITY);
        addSpecies(speciesGhost);

        PotentialMasterMonatomic potentialMaster = new PotentialMasterMonatomic(this); //List(this, 2.0);

        int N = space.D() == 3 ? 256 : 100;  //number of atoms

        double sigma = 1.0;
        double lambda = 1.5;

        //controller and integrator
        box = this.makeBox();
        integrator = new IntegratorHard(this, potentialMaster, box);
        integrator.setTimeStep(0.2);
        integrator.setTemperature(1.0);
        integrator.setIsothermal(false);
        integrator.setThermostat(ThermostatType.ANDERSEN_SCALING);
        integrator.setThermostatNoDrift(true);
        integrator.setThermostatInterval(1);
        P1HardPeriodic nullPotential = new P1HardPeriodic(space, sigma * lambda);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);

        //potentials
        integrator.setNullPotential(nullPotential, species.getLeafType());

        //instantiate several potentials for selection in combo-box
        P2SquareWell potentialSW = new P2SquareWell(space, sigma, lambda, 1.0, true);
        potentialWrapper = new P2HardWrapper(space, potentialSW);
        potentialMaster.addPotential(potentialWrapper, new AtomType[]{species.getLeafType(), species.getLeafType()});

        potentialGhost = new P2DoubleWell(space, 1.0, lambda, 0.0, 0.0);
        potentialMaster.addPotential(potentialGhost, new AtomType[]{species.getLeafType(), speciesGhost.getLeafType()});

        //construct box
        Vector dim = space.makeVector();
        dim.E(space.D() == 3 ? 8 : 13.5);
        box.getBoundary().setBoxSize(dim);
        box.setNMolecules(species, N);
        new ConfigurationLattice(space.D() == 3 ? new LatticeCubicFcc(space) : new LatticeOrthorhombicHexagonal(space), space).initializeCoordinates(box);
        box.setNMolecules(speciesGhost, 1);

        integrator.getEventManager().addListener(new IntegratorListenerAction(new BoxImposePbc(box, space)));
    }
    
    public static void main(String[] args) {
        Space space = Space3D.getInstance();
        if(args.length != 0) {
            try {
                int D = Integer.parseInt(args[0]);
                if (D == 3) {
                    space = Space3D.getInstance();
                }
            } catch(NumberFormatException e) {}
        }
            
        Insertion sim = new Insertion(space);
        sim.getController().actionPerformed();
    }
}
