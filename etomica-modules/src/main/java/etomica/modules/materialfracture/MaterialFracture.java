/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.materialfracture;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorMD;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.BasisOrthorhombicHexagonal;
import etomica.lattice.crystal.PrimitiveGeneral;
import etomica.integrator.IntegratorListenerAction;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncatedForceShifted;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space.Vector;
import etomica.space2d.Space2D;
import etomica.species.SpeciesSpheresMono;

public class MaterialFracture extends Simulation {

    public final Box box;
    public final ConfigurationLattice config;
    public final SpeciesSpheresMono species;
    public final IntegratorVelocityVerlet integrator;
    public final PotentialCalculationForceStress pc;
    public final P2SoftSphericalTruncatedForceShifted pt;
    public final P2LennardJones p2LJ;
    public final P1Tension p1Tension;

    public MaterialFracture() {
        super(Space2D.getInstance());
        PotentialMaster potentialMaster = new PotentialMaster();
        box = new Box(space);
        box.setBoundary(new BoundaryRectangularSlit(0, space));
        box.getBoundary().setBoxSize(space.makeVector(new double[]{90,30}));
        addBox(box);
        integrator = new IntegratorVelocityVerlet(this, potentialMaster, space);
        integrator.setIsothermal(true);
        integrator.setTemperature(300.0);
        integrator.setTimeStep(0.007);
        integrator.setThermostatInterval(100);
        integrator.setThermostat(IntegratorMD.ThermostatType.ANDERSEN);
        integrator.setThermostatNoDrift(true);
        integrator.setBox(box);
        pc = new PotentialCalculationForceStress(space);
        integrator.setForceSum(pc);
        getController().addAction(new ActivityIntegrate(integrator));
        p2LJ = new P2LennardJones(space, 3, 2000);
        pt = new P2SoftSphericalTruncatedForceShifted(space, p2LJ, 7);

        p1Tension = new P1Tension(space); 
        species = new SpeciesSpheresMono(this, space);
        species.setIsDynamic(true);
        ((ElementSimple)species.getLeafType().getElement()).setMass(40);
        addSpecies(species);
        box.setNMolecules(species, 198);

        potentialMaster.addPotential(pt, new AtomType[]{species.getLeafType(), species.getLeafType()});
        potentialMaster.addPotential(p1Tension, new AtomType[]{species.getLeafType()});

        PrimitiveGeneral primitive = new PrimitiveGeneral(space, new Vector[]{space.makeVector(new double[]{Math.sqrt(3),0}), space.makeVector(new double[]{0,1})});
        config = new ConfigurationLattice(new BravaisLatticeCrystal(primitive, new BasisOrthorhombicHexagonal()), space) {
            public void initializeCoordinates(Box aBox) {
                Vector d = space.makeVector();
                d.E(aBox.getBoundary().getBoxSize());
                d.setX(0, 64.7);
                aBox.getBoundary().setBoxSize(d);
                super.initializeCoordinates(aBox);
                d.setX(0, 90);
                aBox.getBoundary().setBoxSize(d);
            }
        };
        config.initializeCoordinates(box);

        integrator.getEventManager().addListener(new IntegratorListenerAction(new BoxImposePbc(box, space)));
    }
}
