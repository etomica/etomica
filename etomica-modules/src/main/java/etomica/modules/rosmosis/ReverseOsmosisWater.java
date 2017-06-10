/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.rosmosis;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.Chlorine;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Sodium;
import etomica.integrator.IntegratorRigidIterative;
import etomica.listener.IntegratorListenerAction;
import etomica.models.water.OrientationCalcWater3P;
import etomica.models.water.P2WaterSPCSoft;
import etomica.models.water.SpeciesWater3POriented;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Dalton;
import etomica.units.Electron;
import etomica.units.Kelvin;

/**
 * Reverse osmosis simulation, based on simulation code by Sohail Murad.
 * 
 * @author Andrew Schultz
 */
public class ReverseOsmosisWater extends Simulation {

    private static final long serialVersionUID = 1L;
    public SpeciesWater3POriented speciesSolvent;
    public SpeciesSpheresMono speciesSodium, speciesChlorine, speciesMembrane;
    public Box box;
    public IntegratorRigidIterative integrator;
    public P2WaterSPCSoft potentialWater;
    public P2LennardJones potentialLJNaNa, potentialLJNaCl, potentialLJClCl;
    public P2Electrostatic potentialQNaNa, potentialQNaCl, potentialQClCl;
    public P2LennardJones potentialLJOCl, potentialLJONa;
    public P2Electrostatic potentialQOCl, potentialQONa;
    public P2Electrostatic potentialQHNa, potentialQHCl;
    public P2LennardJones potentialMM, potentialMO, potentialMCl, potentialMNa;
    public ActivityIntegrate activityIntegrate;
    public ConfigurationMembraneWater configMembrane;
    public P1Tether potentialTether;
    public PotentialCalculationTorqueSumWallForce torqueSum;
    
    public ReverseOsmosisWater(Space space) {
        super(space);
        PotentialMaster potentialMaster = new PotentialMaster(); //List(this, 2.0);
        
        //controller and integrator
	    integrator = new IntegratorRigidIterative(this, potentialMaster, 0.01, Kelvin.UNIT.toSim(298), space);
	    integrator.setIsothermal(true);
        integrator.setThermostatInterval(100);
        integrator.setTimeStep(0.004);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);

        //solute (1)
        speciesSodium = new SpeciesSpheresMono(space, Sodium.INSTANCE);
        speciesSodium.setIsDynamic(true);
        addSpecies(speciesSodium);
        
        speciesChlorine = new SpeciesSpheresMono(space, Chlorine.INSTANCE);
        speciesChlorine.setIsDynamic(true);
        addSpecies(speciesChlorine);
        
        //solvent (2)
        speciesSolvent = new SpeciesWater3POriented(space, true);
        addSpecies(speciesSolvent);
        integrator.setOrientationCalc(speciesSolvent, new OrientationCalcWater3P(space));

        //membrane
        speciesMembrane = new SpeciesSpheresMono(this, space);
        speciesMembrane.setIsDynamic(true);
        ((ElementSimple)speciesMembrane.getLeafType().getElement()).setMass(Dalton.UNIT.toSim(80));
        addSpecies(speciesMembrane);
        
        /*
         * Sodium and chloride potential parameters from 
         * I G. Tironi et al., J. Chem. Phys. 102 (1995), pp. 5451
         */
        
        double epsSodium = 6.177425;
        double epsChlorine = 44.5586;
        double sigSodium = 2.57155;
        double sigChlorine = 4.44804;
        double chargeSodium = Electron.UNIT.toSim(1);
        double chargeChlorine = Electron.UNIT.toSim(-1);
        double chargeOxygen = Electron.UNIT.toSim(-0.82);
        double chargeHydrogen = Electron.UNIT.toSim(0.41);
        double epsMembrane = Kelvin.UNIT.toSim(12.5);
        double sigMembrane = 3;
        
        double xSize = 44;// 80 originally
        double yzSize = 18;       // 28 originally
        double rCut = 0.49*yzSize;
        double switchFac = 0.7;

        AtomType oType = speciesSolvent.getOxygenType();
        AtomType hType = speciesSolvent.getHydrogenType();
        AtomType naType = speciesSodium.getLeafType();
        AtomType clType = speciesChlorine.getLeafType();
        AtomType mType = speciesMembrane.getLeafType();
        
        potentialWater = new P2WaterSPCSoft(space);
        P2MoleculeSoftTruncatedSwitched pWaterTrunc = new P2MoleculeSoftTruncatedSwitched(potentialWater, rCut, space);
        pWaterTrunc.setSwitchFac(switchFac);
        potentialMaster.addPotential(pWaterTrunc, new ISpecies[]{speciesSolvent, speciesSolvent});
        double epsOxygen = potentialWater.getEpsilon();
        double sigOxygen = potentialWater.getSigma();
        
        potentialLJNaNa = new P2LennardJones(space, sigSodium, epsSodium);
	    P2SoftSphericalTruncatedSwitched pTrunc = new P2SoftSphericalTruncatedSwitched(space, potentialLJNaNa, rCut);
	    pTrunc.setSwitchFac(switchFac);
        potentialMaster.addPotential(pTrunc, new AtomType[]{naType, naType});

        potentialLJClCl = new P2LennardJones(space, sigChlorine, epsChlorine);
        pTrunc = new P2SoftSphericalTruncatedSwitched(space, potentialLJClCl, rCut);
        pTrunc.setSwitchFac(switchFac);
        potentialMaster.addPotential(pTrunc, new AtomType[]{clType, clType});

        potentialLJNaCl = new P2LennardJones(space, 0.5*(sigChlorine+sigSodium), Math.sqrt(epsChlorine*epsSodium));
        pTrunc = new P2SoftSphericalTruncatedSwitched(space, potentialLJNaCl, rCut);
        pTrunc.setSwitchFac(switchFac);
        potentialMaster.addPotential(pTrunc, new AtomType[]{naType, clType});

        potentialQNaNa = new P2Electrostatic(space);
        potentialQNaNa.setCharge1(chargeSodium);
        potentialQNaNa.setCharge2(chargeSodium);
        pTrunc = new P2SoftSphericalTruncatedSwitched(space, potentialQNaNa, rCut);
        pTrunc.setSwitchFac(switchFac);
        potentialMaster.addPotential(pTrunc, new AtomType[]{naType, naType});

        potentialQClCl = new P2Electrostatic(space);
        potentialQClCl.setCharge1(chargeChlorine);
        potentialQClCl.setCharge2(chargeChlorine);
        pTrunc = new P2SoftSphericalTruncatedSwitched(space, potentialQClCl, rCut);
        pTrunc.setSwitchFac(switchFac);
        potentialMaster.addPotential(pTrunc, new AtomType[]{clType, clType});

        potentialQNaCl = new P2Electrostatic(space);
        potentialQNaCl.setCharge1(chargeSodium);
        potentialQNaCl.setCharge2(chargeChlorine);
        pTrunc = new P2SoftSphericalTruncatedSwitched(space, potentialQNaCl, rCut);
        pTrunc.setSwitchFac(switchFac);
        potentialMaster.addPotential(pTrunc, new AtomType[]{naType, clType});

        potentialLJOCl = new P2LennardJones(space, 0.5*(sigOxygen+sigChlorine), Math.sqrt(epsOxygen*epsChlorine));
        pTrunc = new P2SoftSphericalTruncatedSwitched(space, potentialLJOCl, rCut);
        pTrunc.setSwitchFac(switchFac);
        potentialMaster.addPotential(pTrunc, new AtomType[]{clType, oType});

        potentialLJONa = new P2LennardJones(space, 0.5*(sigOxygen+sigSodium), Math.sqrt(epsOxygen*epsSodium));
        pTrunc = new P2SoftSphericalTruncatedSwitched(space, potentialLJONa, rCut);
        pTrunc.setSwitchFac(switchFac);
        potentialMaster.addPotential(pTrunc, new AtomType[]{naType, oType});

        potentialQONa = new P2Electrostatic(space);
        potentialQONa.setCharge1(chargeOxygen);
        potentialQONa.setCharge2(chargeSodium);
        pTrunc = new P2SoftSphericalTruncatedSwitched(space, potentialQONa, rCut);
        pTrunc.setSwitchFac(switchFac);
        potentialMaster.addPotential(pTrunc, new AtomType[]{oType, naType});
        
        potentialQOCl = new P2Electrostatic(space);
        potentialQOCl.setCharge1(chargeOxygen);
        potentialQOCl.setCharge2(chargeChlorine);
        pTrunc = new P2SoftSphericalTruncatedSwitched(space, potentialQOCl, rCut);
        pTrunc.setSwitchFac(switchFac);
        potentialMaster.addPotential(pTrunc, new AtomType[]{oType, clType});

        potentialQHNa = new P2Electrostatic(space);
        potentialQHNa.setCharge1(chargeHydrogen);
        potentialQHNa.setCharge2(chargeSodium);
        pTrunc = new P2SoftSphericalTruncatedSwitched(space, potentialQHNa, rCut);
        pTrunc.setSwitchFac(switchFac);
        potentialMaster.addPotential(pTrunc, new AtomType[]{hType, naType});
        
        potentialQHCl = new P2Electrostatic(space);
        potentialQHCl.setCharge1(chargeHydrogen);
        potentialQHCl.setCharge2(chargeChlorine);
        pTrunc = new P2SoftSphericalTruncatedSwitched(space, potentialQHCl, rCut);
        pTrunc.setSwitchFac(switchFac);
        potentialMaster.addPotential(pTrunc, new AtomType[]{hType, clType});

        potentialMM = new P2LennardJones(space, sigMembrane, epsMembrane);
        pTrunc = new P2SoftSphericalTruncatedSwitched(space, potentialMM, rCut);
        pTrunc.setSwitchFac(switchFac);
        potentialMaster.addPotential(pTrunc, new AtomType[]{mType, mType});
        
        potentialMO = new P2LennardJones(space, 0.5*(sigMembrane+sigOxygen), Math.sqrt(epsMembrane*epsOxygen));
        pTrunc = new P2SoftSphericalTruncatedSwitched(space, potentialMO, rCut);
        pTrunc.setSwitchFac(switchFac);
        potentialMaster.addPotential(pTrunc, new AtomType[]{mType, oType});
        
        potentialMNa = new P2LennardJones(space, 0.5*(sigMembrane+sigSodium), Math.sqrt(epsMembrane*epsSodium));
        pTrunc = new P2SoftSphericalTruncatedSwitched(space, potentialMNa, rCut);
        pTrunc.setSwitchFac(switchFac);
        potentialMaster.addPotential(pTrunc, new AtomType[]{mType, naType});

        potentialMCl = new P2LennardJones(space, 0.5*(sigMembrane+sigChlorine), Math.sqrt(epsMembrane*epsChlorine));
        pTrunc = new P2SoftSphericalTruncatedSwitched(space, potentialMCl, rCut);
        pTrunc.setSwitchFac(switchFac);
        potentialMaster.addPotential(pTrunc, new AtomType[]{mType, clType});

        //construct box
	    box = new Box(new BoundaryRectangularPeriodic(space, 15), space);
        addBox(box);
        Vector dim = space.makeVector();
        dim.E(new double[]{xSize, yzSize, yzSize});
        box.getBoundary().setBoxSize(dim);
        configMembrane = new ConfigurationMembraneWater(this, space);
        configMembrane.setMembraneDim(0);
        configMembrane.setMembraneThickness(2*4.0);
        configMembrane.setNumMembraneLayers(2);
        configMembrane.setMembraneWidth(2);
        double density = 0.01;
        configMembrane.setSolutionChamberDensity(density);
        configMembrane.setSolventChamberDensity(density);
        configMembrane.setSoluteMoleFraction(0.5);
        configMembrane.setSpeciesMembrane(speciesMembrane);
        configMembrane.setSpeciesSolute1(speciesSodium);
        configMembrane.setSpeciesSolute2(speciesChlorine);
        configMembrane.setSpeciesSolvent(speciesSolvent);
        configMembrane.initializeCoordinates(box);

        potentialTether = new P1Tether(box, speciesMembrane, space);
        potentialTether.setEpsilon(20000*298/125);
        potentialMaster.addPotential(potentialTether, new AtomType[]{speciesMembrane.getLeafType()});
        
        integrator.setBox(box);

        BoxImposePbc pbc = new BoxImposePbc(box, space);
        pbc.setApplyToMolecules(true);
        integrator.getEventManager().addListener(new IntegratorListenerAction(pbc));
        
        torqueSum = new PotentialCalculationTorqueSumWallForce(potentialTether);
        integrator.setTorqueSum(torqueSum);
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
            
        ReverseOsmosisWater sim = new ReverseOsmosisWater(space);
        sim.getController().actionPerformed();
    }//end of main
    
}
