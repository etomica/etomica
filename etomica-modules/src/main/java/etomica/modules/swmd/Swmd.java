/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.swmd;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMD.ThermostatType;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.integrator.IntegratorListenerAction;
import etomica.potential.P1HardPeriodic;
import etomica.potential.P2HardWrapper;
import etomica.potential.P2SquareWell;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.*;

public class Swmd extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public SpeciesSpheresMono species;
    public Box box;
    public IntegratorHard integrator;
    public P2HardWrapper potentialWrapper;
    public ActivityIntegrate activityIntegrate;
    
    public Swmd(Space _space) {
        super(_space);
        PotentialMasterMonatomic potentialMaster = new PotentialMasterMonatomic(this); //List(this, 2.0);

        int N = space.D() == 3 ? 256 : 100;  //number of atoms

        double sigma = 4.0;
        double lambda = 2.0;

        //controller and integrator
        box = new Box(space);
        integrator = new IntegratorHard(this, potentialMaster, space, box);
        integrator.setTimeStep(1.0);
        integrator.setTemperature(Kelvin.UNIT.toSim(300));
        integrator.setIsothermal(false);
        integrator.setThermostat(ThermostatType.ANDERSEN_SINGLE);
        integrator.setThermostatInterval(1);
        P1HardPeriodic nullPotential = new P1HardPeriodic(space, sigma * lambda);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);

        //species and potentials
        species = new SpeciesSpheresMono(this, space);//index 1
        species.setIsDynamic(true);
        ((ElementSimple) species.getLeafType().getElement()).setMass(Dalton.UNIT.toSim(space.D() == 3 ? 131 : 40));
        addSpecies(species);
        integrator.setNullPotential(nullPotential, species.getLeafType());

        //instantiate several potentials for selection in combo-box
        P2SquareWell potentialSW = new P2SquareWell(space, sigma, lambda, new UnitRatio(Joule.UNIT, Mole.UNIT).toSim(space.D() == 3 ? 1000 : 1500), true);
        potentialWrapper = new P2HardWrapper(space, potentialSW);
        potentialMaster.addPotential(potentialWrapper, new AtomType[]{species.getLeafType(), species.getLeafType()});

        //construct box
        addBox(box);
        Vector dim = space.makeVector();
        dim.E(space.D() == 3 ? 30 : 50);
        box.getBoundary().setBoxSize(dim);
        box.setNMolecules(species, N);
        new ConfigurationLattice(space.D() == 3 ? new LatticeCubicFcc(space) : new LatticeOrthorhombicHexagonal(space), space).initializeCoordinates(box);
        integrator.setBox(box);

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
            
        Swmd sim = new Swmd(space);
        sim.getController().actionPerformed();
    }
}
