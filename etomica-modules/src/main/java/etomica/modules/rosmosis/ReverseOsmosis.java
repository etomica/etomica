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
import etomica.integrator.IntegratorMD.ThermostatType;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncatedShifted;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.Space;
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
public class ReverseOsmosis extends Simulation {

    private static final long serialVersionUID = 1L;
    public SpeciesGeneral speciesSolvent, speciesSolute, speciesMembrane;
    public Box box;
    public IntegratorVelocityVerlet integrator;
    public P2LennardJones potential11, potential12, potential22;
    public P2LennardJones potentialMM, potentialM1, potentialM2;

    public ConfigurationMembrane configMembrane;
    public P1Tether potentialTether;
    public PotentialCalculationForceSumWallForce forceSum;
    
    public ReverseOsmosis() {
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

        PotentialMaster potentialMaster = new PotentialMasterMonatomic(this); //List(this, 2.0);

        //controller and integrator
        box = this.makeBox();
        integrator = new IntegratorVelocityVerlet(potentialMaster, getRandom(), 0.01, Kelvin.UNIT.toSim(125), box);
        integrator.setIsothermal(false);
        integrator.setThermostat(ThermostatType.ANDERSEN_SINGLE);
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
        potential11 = new P2LennardJones(space, sigSolute, epsSolute);
        P2SoftSphericalTruncatedShifted pTrunc = new P2SoftSphericalTruncatedShifted(space, potential11, rCut);
        potentialMaster.addPotential(pTrunc, new AtomType[]{speciesSolute.getLeafType(), speciesSolute.getLeafType()});

        potential22 = new P2LennardJones(space, sigSolvent, epsSolvent);
        pTrunc = new P2SoftSphericalTruncatedShifted(space, potential22, rCut);
        potentialMaster.addPotential(pTrunc, new AtomType[]{speciesSolvent.getLeafType(), speciesSolvent.getLeafType()});

        potential12 = new P2LennardJones(space, 0.5 * (sigSolvent + sigSolute), Math.sqrt(epsSolvent * epsSolute));
        pTrunc = new P2SoftSphericalTruncatedShifted(space, potential12, rCut);
        potentialMaster.addPotential(pTrunc, new AtomType[]{speciesSolvent.getLeafType(), speciesSolute.getLeafType()});

        potentialMM = new P2LennardJones(space, sigMembrane, epsMembrane);
        pTrunc = new P2SoftSphericalTruncatedShifted(space, potentialMM, rCut);
        potentialMaster.addPotential(pTrunc, new AtomType[]{speciesMembrane.getLeafType(), speciesMembrane.getLeafType()});

        potentialM1 = new P2LennardJones(space, 0.5 * (sigMembrane + sigSolute), Math.sqrt(epsMembrane * epsSolute));
        pTrunc = new P2SoftSphericalTruncatedShifted(space, potentialM1, rCut);
        potentialMaster.addPotential(pTrunc, new AtomType[]{speciesMembrane.getLeafType(), speciesSolute.getLeafType()});

        potentialM2 = new P2LennardJones(space, 0.5 * (sigMembrane + sigSolvent), Math.sqrt(epsMembrane * epsSolvent));
        pTrunc = new P2SoftSphericalTruncatedShifted(space, potentialM2, rCut);
        potentialMaster.addPotential(pTrunc, new AtomType[]{speciesMembrane.getLeafType(), speciesSolvent.getLeafType()});


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
        potentialMaster.addPotential(potentialTether, new AtomType[]{speciesMembrane.getLeafType()});

//        integrator.addIntervalAction(potentialMaster.getNeighborManager(box));
//        integrator.addNonintervalListener(potentialMaster.getNeighborManager(box));
        integrator.getEventManager().addListener(new IntegratorListenerAction(new BoxImposePbc(box, space)));
        forceSum = new PotentialCalculationForceSumWallForce(potentialTether, integrator.getForces());
        integrator.setForceSum(forceSum);
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
            
        ReverseOsmosis sim = new ReverseOsmosis();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, Long.MAX_VALUE));
    }//end of main
    
}
