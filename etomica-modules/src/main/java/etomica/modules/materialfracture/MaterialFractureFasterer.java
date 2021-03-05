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
import etomica.integrator.IntegratorMDFasterer;
import etomica.integrator.IntegratorVelocityVerletFasterer;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.BasisOrthorhombicHexagonal;
import etomica.lattice.crystal.PrimitiveGeneral;
import etomica.potential.BondingInfo;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncatedForceShifted;
import etomica.potential.PotentialMasterFasterer;
import etomica.potential.compute.PotentialCallback;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputeField;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space.Vector;
import etomica.space2d.Space2D;
import etomica.species.SpeciesGeneral;

public class MaterialFractureFasterer extends Simulation {

    public final Box box;
    public final ConfigurationLattice config;
    public final SpeciesGeneral species;
    public final IntegratorVelocityVerletFasterer integrator;
    public final P2SoftSphericalTruncatedForceShifted pt;
    public final P2LennardJones p2LJ;
    public final P1Tension p1Tension;
    public double wallForce;

    public MaterialFractureFasterer() {
        super(Space2D.getInstance());

        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        ((ElementSimple) species.getLeafType().getElement()).setMass(40);
        addSpecies(species);

        box = this.makeBox(new BoundaryRectangularSlit(0, space));
        PotentialMasterFasterer potentialMaster = new PotentialMasterFasterer(getSpeciesManager(), box, BondingInfo.noBonding());
        PotentialComputeField pcField = new PotentialComputeField(getSpeciesManager(), box) {
            public double computeAll(boolean doForces, PotentialCallback pc) {
                double u = super.computeAll(doForces, pc);
                wallForce = 0;
                for (Vector f : forces) {
                    wallForce += Math.abs(f.getX(0));
                }
                return u;
            }
        };
        PotentialComputeAggregate pcAgg = new PotentialComputeAggregate(potentialMaster, pcField);

        box.getBoundary().setBoxSize(Vector.of(new double[]{90, 30}));
        integrator = new IntegratorVelocityVerletFasterer(pcAgg, this.getRandom(), 0.05, 1.0, box);
        integrator.setIsothermal(true);
        integrator.setTemperature(300.0);
        integrator.setTimeStep(0.007);
        integrator.setThermostatInterval(100);
        integrator.setThermostat(IntegratorMDFasterer.ThermostatType.ANDERSEN);
        integrator.setThermostatNoDrift(true);
        getController().addActivity(new ActivityIntegrate(integrator));
        p2LJ = new P2LennardJones(space, 3, 2000);
        pt = new P2SoftSphericalTruncatedForceShifted(space, p2LJ, 7);

        box.setNMolecules(species, 198);

        potentialMaster.setPairPotential(species.getLeafType(), species.getLeafType(), pt);
        p1Tension = new P1Tension(space);
        pcField.setFieldPotential(species.getLeafType(), p1Tension);

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
