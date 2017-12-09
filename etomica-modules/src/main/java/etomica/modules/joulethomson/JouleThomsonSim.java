/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.joulethomson;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.graphics.DisplayBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorGear4NPH;
import etomica.integrator.IntegratorMD;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.lattice.SpaceLattice;
import etomica.listener.IntegratorListenerAction;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.*;

public class JouleThomsonSim extends Simulation {
    
    private static final long serialVersionUID = 1L;
    IntegratorGear4NPH integrator;
    IntegratorMD integratorNVE;
    SpeciesSpheresMono species;
    P2LennardJones potential;
    Box box;
    IntegratorJT integratorJT;
    ActivityIntegrate activityIntegrate;
    Configuration config;

    public JouleThomsonSim() {this(Space2D.getInstance());}
    public JouleThomsonSim(Space space) {
        super(space);
        PotentialMaster potentialMaster = new PotentialMasterMonatomic(this);
        int nAtoms = (space.D() < 3) ? 50 : 64;
        double sigma = 3.0;
        
        //integrator
        integratorNVE = new IntegratorVelocityVerlet(this, potentialMaster, space);
        integrator = new IntegratorGear4NPH(this, potentialMaster, space);
        integrator.setRelaxationRateP(500.);
        integrator.setRelaxationRateH(300.);
        integratorNVE.setTemperature(Kelvin.UNIT.toSim(300));
        integrator.setTemperature(Kelvin.UNIT.toSim(300));
        final Unit pUnit;
        if (space.D() == 2) {
            Unit[] units = new Unit[] {Bar.UNIT, new PrefixedUnit(Prefix.NANO, Meter.UNIT)};
            double[] exponents = new double[] {1.0, 1.0};
            pUnit = new CompoundUnit(units, exponents);
        }
        else {
            pUnit = Bar.UNIT;
        }
        integrator.setTargetP(pUnit.toSim(100.0));
	    
	    //species and potential
	    species = new SpeciesSpheresMono(this, space);
	    species.setIsDynamic(true);
	    
        ((ElementSimple)species.getLeafType().getElement()).setMass(40);
        addSpecies(species);
	    potential = new P2LennardJones(space, sigma, Kelvin.UNIT.toSim(300));
        potentialMaster.addPotential(potential, new AtomType[]{species.getLeafType(), species.getLeafType()});
        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, nAtoms);
        
        SpaceLattice lattice;
        if (space.D() == 2) {
            lattice = new LatticeOrthorhombicHexagonal(space);
        }
        else {
            lattice = new LatticeCubicFcc(space);
        }
        config = new ConfigurationLattice(lattice, space);
        config.initializeCoordinates(box);
        
        integratorJT = new IntegratorJT(getRandom(), integrator, integratorNVE);
        integratorJT.setTemperature(Kelvin.UNIT.toSim(300));
        integratorJT.getEventManager().addListener(new IntegratorListenerAction(new BoxImposePbc(box, space)));
        integrator.setBox(box);
        integratorNVE.setBox(box);

        integrator.setTimeStep(0.001);
        integratorNVE.setTimeStep(0.005);

        activityIntegrate = new ActivityIntegrate(integratorJT);
        getController().addAction(activityIntegrate);
    }
    
    public static void main(String[] args) {
    	Space sp = Space3D.getInstance();
        JouleThomsonSim sim = new JouleThomsonSim(sp);
        SimulationGraphic simGraphic = new SimulationGraphic(sim);
        simGraphic.makeAndDisplayFrame();
        final DisplayBox displayBox = simGraphic.getDisplayBox(sim.box);
        sim.activityIntegrate.setSleepPeriod(10);

        sim.integratorJT.getEventManager().addListener(new IntegratorListenerAction(simGraphic.getPaintAction(sim.box)));

    }
}
