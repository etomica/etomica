/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.rosmosis;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMDFasterer;
import etomica.integrator.IntegratorVelocityVerletFasterer;
import etomica.potential.*;
import etomica.potential.compute.PotentialCallback;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputeField;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.units.Dalton;
import etomica.units.Kelvin;

/**
 * Reverse osmosis simulation, based on simulation code by Sohail Murad.
 *
 * @author Andrew Schultz
 */
public class ReverseOsmosisFasterer extends Simulation implements MeterOsmoticPressureFasterer.WallForceSource {

    public SpeciesGeneral speciesSolvent, speciesSolute, speciesMembrane;
    public Box box;
    public IntegratorVelocityVerletFasterer integrator;
    public P2LennardJones potential11, potential12, potential22;
    public P2LennardJones potentialMM, potentialM1, potentialM2;

    public ConfigurationMembrane configMembrane;
    public P1Tether potentialTether;
    public double wallForce;

    public ReverseOsmosisFasterer() {
        super(Space3D.getInstance());

        //solute (1)
        speciesSolute = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        ((ElementSimple) speciesSolute.getLeafType().getElement()).setMass(Dalton.UNIT.toSim(40));
        addSpecies(speciesSolute);

        //solvent (2)
        speciesSolvent = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        ((ElementSimple) speciesSolvent.getLeafType().getElement()).setMass(Dalton.UNIT.toSim(40));
        addSpecies(speciesSolvent);

        //membrane
        speciesMembrane = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        ((ElementSimple) speciesMembrane.getLeafType().getElement()).setMass(Dalton.UNIT.toSim(80));
        addSpecies(speciesMembrane);
        box = this.makeBox();

        PotentialMasterFasterer potentialMaster = new PotentialMasterFasterer(getSpeciesManager(), box, BondingInfo.noBonding());
        PotentialComputeField pcField = new PotentialComputeField(getSpeciesManager(), box) {
            @Override
            public double computeAll(boolean doForces, PotentialCallback pc) {
                double u = super.computeAll(doForces, pc);
                if (doForces) {
                    wallForce = 0;
                    int n = box.getLeafList().size();
                    for (int i = 0; i < n; i++) {
                        double fi = forces[i].getX(0);
                        if (box.getLeafList().get(i).getPosition().getX(0) > 0) fi = -fi;
                        wallForce += fi;
                    }
                }
                return u;
            }
        };
        PotentialComputeAggregate pcAgg = new PotentialComputeAggregate(potentialMaster, pcField);

        //controller and integrator
        integrator = new IntegratorVelocityVerletFasterer(pcAgg, getRandom(), 0.01, Kelvin.UNIT.toSim(125), box);
        integrator.setIsothermal(false);
        integrator.setThermostat(IntegratorMDFasterer.ThermostatType.ANDERSEN_SINGLE);
        integrator.setThermostatInterval(1);
        integrator.setTimeStep(0.02);

        getController().addActivity(new ActivityIntegrate(integrator));

        double epsSolute = Kelvin.UNIT.toSim(125.0);
        double sigSolute = 3.5;
        double epsSolvent = Kelvin.UNIT.toSim(125.0);
        double sigSolvent = 0.5 * 3.5;
        double epsMembrane = Kelvin.UNIT.toSim(12.5);
        double sigMembrane = 0.988 * 3.5; // ???

        double xSize = 66 + 2.0 / 3.0;// 80 originally
        double yzSize = 21;       // 28 originally
        double rCut = 0.5 * yzSize;

        //instantiate several potentials for selection in combo-box
        TruncationFactory tf = new TruncationFactoryShift(space, rCut);
        potential11 = new P2LennardJones(space, sigSolute, epsSolute);
        potentialMaster.setPairPotential(speciesSolute.getLeafType(), speciesSolute.getLeafType(), tf.make(potential11));

        potential22 = new P2LennardJones(space, sigSolvent, epsSolvent);
        potentialMaster.setPairPotential(speciesSolvent.getLeafType(), speciesSolvent.getLeafType(), tf.make(potential22));

        potential12 = new P2LennardJones(space, 0.5 * (sigSolvent + sigSolute), Math.sqrt(epsSolvent * epsSolute));
        potentialMaster.setPairPotential(speciesSolvent.getLeafType(), speciesSolute.getLeafType(), tf.make(potential12));

        potentialMM = new P2LennardJones(space, sigMembrane, epsMembrane);
        potentialMaster.setPairPotential(speciesMembrane.getLeafType(), speciesMembrane.getLeafType(), potentialMM);

        potentialM1 = new P2LennardJones(space, 0.5 * (sigMembrane + sigSolute), Math.sqrt(epsMembrane * epsSolute));
        potentialMaster.setPairPotential(speciesMembrane.getLeafType(), speciesSolute.getLeafType(), tf.make(potentialM1));

        potentialM2 = new P2LennardJones(space, 0.5 * (sigMembrane + sigSolvent), Math.sqrt(epsMembrane * epsSolvent));
        potentialMaster.setPairPotential(speciesMembrane.getLeafType(), speciesSolvent.getLeafType(), tf.make(potentialM2));


        //construct box
        Vector dim = space.makeVector();
        dim.E(new double[]{xSize, yzSize, yzSize});
        box.getBoundary().setBoxSize(dim);
        configMembrane = new ConfigurationMembrane(this, space);
        configMembrane.setMembraneDim(0);
        configMembrane.setMembraneThicknessPerLayer(10.0 / 3.0);
        configMembrane.setNumMembraneLayers(2);
        configMembrane.setMembraneWidth(3);
        double density = 0.525 / Math.pow(sigSolute, 3);
        configMembrane.setSolutionChamberDensity(density);
        configMembrane.setSolventChamberDensity(density);
        configMembrane.setSpeciesMembrane(speciesMembrane);
        configMembrane.setSpeciesSolute(speciesSolute);
        configMembrane.setSpeciesSolvent(speciesSolvent);
        configMembrane.initializeCoordinates(box);

        potentialTether = new P1Tether(box, speciesMembrane, space);
        potentialTether.setEpsilon(20000);
        pcField.setFieldPotential(speciesMembrane.getLeafType(), potentialTether);

        integrator.getEventManager().addListener(new IntegratorListenerAction(new BoxImposePbc(box, space)));
    }

    public double getWallForce() {
        return wallForce;
    }

    public static void main(String[] args) {
        ReverseOsmosisFasterer sim = new ReverseOsmosisFasterer();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, Long.MAX_VALUE));
    }//end of main

}
