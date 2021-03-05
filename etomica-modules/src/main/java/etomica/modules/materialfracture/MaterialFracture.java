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
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMD;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.BasisOrthorhombicHexagonal;
import etomica.lattice.crystal.PrimitiveGeneral;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncatedForceShifted;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space.Vector;
import etomica.space2d.Space2D;
import etomica.species.SpeciesGeneral;

public class MaterialFracture extends Simulation {

    public final Box box;
    public final ConfigurationLattice config;
    public final SpeciesGeneral species;
    public final IntegratorVelocityVerlet integrator;
    public final PotentialCalculationForceStress pc;
    public final P2SoftSphericalTruncatedForceShifted pt;
    public final P2LennardJones p2LJ;
    public final P1Tension p1Tension;

    public MaterialFracture() {
        super(Space2D.getInstance());

        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        ((ElementSimple) species.getLeafType().getElement()).setMass(40);
        addSpecies(species);

        PotentialMaster potentialMaster = new PotentialMaster();
        box = this.makeBox(new BoundaryRectangularSlit(0, space));
        box.getBoundary().setBoxSize(Vector.of(new double[]{90, 30}));
        integrator = new IntegratorVelocityVerlet(potentialMaster, random, 0.05, 1.0, box);
        integrator.setIsothermal(true);
        integrator.setTemperature(300.0);
        integrator.setTimeStep(0.007);
        integrator.setThermostatInterval(100);
        integrator.setThermostat(IntegratorMD.ThermostatType.ANDERSEN);
        integrator.setThermostatNoDrift(true);
        pc = new PotentialCalculationForceStress(space);
        integrator.setForceSum(pc);
        getController().addActivity(new ActivityIntegrate(integrator));
        p2LJ = new P2LennardJones(space, 3, 2000);
        pt = new P2SoftSphericalTruncatedForceShifted(space, p2LJ, 7);

        p1Tension = new P1Tension(space);
        box.setNMolecules(species, 198);

        potentialMaster.addPotential(pt, new AtomType[]{species.getLeafType(), species.getLeafType()});
        potentialMaster.addPotential(p1Tension, new AtomType[]{species.getLeafType()});

        PrimitiveGeneral primitive = new PrimitiveGeneral(space, new Vector[]{Vector.of(new double[]{Math.sqrt(3), 0}), Vector.of(new double[]{0, 1})});
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
