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
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.integrator.ShakeListener;
import etomica.models.water.ConformationWater3P;
import etomica.models.water.P2WaterSPC;
import etomica.models.water.SpeciesWater3P;
import etomica.potential.*;
import etomica.potential.compute.PotentialCallback;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputeEwaldFourier;
import etomica.potential.compute.PotentialComputeField;
import etomica.potential.ewald.P2Ewald1Real;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesAgentManager;
import etomica.species.SpeciesGeneral;
import etomica.units.Dalton;
import etomica.units.Electron;
import etomica.units.Kelvin;

/**
 * Reverse osmosis simulation, based on simulation code by Sohail Murad.
 *
 * @author Andrew Schultz
 */
public class ReverseOsmosisWater extends Simulation implements MeterOsmoticPressure.WallForceSource {

    public SpeciesGeneral speciesSolvent;
    public SpeciesGeneral speciesSodium, speciesChlorine, speciesMembrane;
    public Box box;
    public IntegratorVelocityVerlet integrator;
    public P2LennardJones potentialLJOO;
    public P2LennardJones potentialLJNaNa, potentialLJNaCl, potentialLJClCl;
    public P2Ewald1Real potentialQNaNa, potentialQNaCl, potentialQClCl;
    public P2LennardJones potentialLJOCl, potentialLJONa;
    public P2Ewald1Real potentialQOCl, potentialQONa;
    public P2Ewald1Real potentialQHNa, potentialQHCl;
    public P2LennardJones potentialMM, potentialMO, potentialMCl, potentialMNa;
    public double wallForce;

    public ConfigurationMembraneWater configMembrane;
    public P1Tether potentialTether;

    public ReverseOsmosisWater() {
        super(Space3D.getInstance());

        //solute (1)
        speciesSodium = SpeciesGeneral.monatomic(space, AtomType.element(Sodium.INSTANCE), true);
        addSpecies(speciesSodium);

        speciesChlorine = SpeciesGeneral.monatomic(space, AtomType.element(Chlorine.INSTANCE), true);
        addSpecies(speciesChlorine);

        //solvent (2)
        speciesSolvent = SpeciesWater3P.create(true);
        addSpecies(speciesSolvent);

        //membrane
        speciesMembrane = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        ((ElementSimple) speciesMembrane.getLeafType().getElement()).setMass(Dalton.UNIT.toSim(80));
        addSpecies(speciesMembrane);

        //construct box
        box = this.makeBox(new BoundaryRectangularPeriodic(space, 15));
        Vector dim = space.makeVector();
        double xSize = 44;// 80 originally
        double yzSize = 18;       // 28 originally
        dim.E(new double[]{xSize, yzSize, yzSize});
        box.getBoundary().setBoxSize(dim);

        configMembrane = new ConfigurationMembraneWater(this, space);
        configMembrane.setMembraneDim(0);
        configMembrane.setMembraneThickness(2 * 4.0);
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

        PotentialMaster potentialMaster = new PotentialMaster(getSpeciesManager(), box, BondingInfo.noBonding()); //List(this, 2.0);

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
        double chargeOxygen = P2WaterSPC.chargeO;
        double chargeHydrogen = P2WaterSPC.chargeH;
        double epsMembrane = Kelvin.UNIT.toSim(12.5);
        double sigMembrane = 3;

        double rCut = 0.49 * yzSize;

        AtomType oType = speciesSolvent.getTypeByName("O");
        AtomType hType = speciesSolvent.getTypeByName("H");
        AtomType naType = speciesSodium.getLeafType();
        AtomType clType = speciesChlorine.getLeafType();
        AtomType mType = speciesMembrane.getLeafType();

        PotentialComputeEwaldFourier ewaldFourier = new PotentialComputeEwaldFourier(getSpeciesManager(), box);
        PotentialComputeEwaldFourier.EwaldParams ewaldParams = ewaldFourier.getOptimalParams(2.5, 0);
        double alpha = ewaldParams.alpha;
        ewaldFourier.setkCut(ewaldParams.kCut);
        ewaldFourier.setAlpha(alpha);
        ewaldFourier.setCharge(oType, chargeOxygen);
        ewaldFourier.setCharge(hType, chargeHydrogen);
        ewaldFourier.setCharge(naType, chargeSodium);
        ewaldFourier.setCharge(clType, chargeChlorine);

        PotentialMasterBonding pmBonding = ewaldFourier.makeIntramolecularCorrection();

        TruncationFactoryForceShift tf = new TruncationFactoryForceShift(space, rCut);
        potentialLJOO = new P2LennardJones(space, P2WaterSPC.sigmaOO, P2WaterSPC.epsilonOO);
        P2Ewald1Real p2OOqq = new P2Ewald1Real(chargeOxygen * chargeOxygen, alpha);
        Potential2Soft p2OO = tf.make(potentialLJOO, p2OOqq);
        Potential2Soft p2OHqq = P2Ewald1Real.makeTruncated(chargeOxygen * chargeHydrogen, ewaldParams.alpha, tf);
        Potential2Soft p2HHqq = P2Ewald1Real.makeTruncated(chargeHydrogen * chargeHydrogen, ewaldParams.alpha, tf);
        potentialMaster.setPairPotential(oType, oType, p2OO);
        potentialMaster.setPairPotential(oType, hType, p2OHqq);
        potentialMaster.setPairPotential(hType, hType, p2HHqq);
        double epsOxygen = P2WaterSPC.epsilonOO;
        double sigOxygen = P2WaterSPC.sigmaOO;

        potentialLJNaNa = new P2LennardJones(space, sigSodium, epsSodium);
        potentialQNaNa = new P2Ewald1Real(chargeSodium * chargeSodium, ewaldParams.alpha);
        potentialMaster.setPairPotential(naType, naType, tf.make(potentialLJNaNa, potentialQNaNa));

        potentialLJClCl = new P2LennardJones(space, sigChlorine, epsChlorine);
        potentialQClCl = new P2Ewald1Real(chargeChlorine * chargeChlorine, ewaldParams.alpha);
        potentialMaster.setPairPotential(clType, clType, tf.make(potentialLJClCl, potentialQClCl));

        potentialLJNaCl = new P2LennardJones(space, 0.5 * (sigChlorine + sigSodium), Math.sqrt(epsChlorine * epsSodium));
        potentialQNaCl = new P2Ewald1Real(chargeSodium * chargeChlorine, ewaldParams.alpha);
        potentialMaster.setPairPotential(naType, clType, tf.make(potentialLJNaCl, potentialQNaCl));

        potentialLJOCl = new P2LennardJones(space, 0.5 * (sigOxygen + sigChlorine), Math.sqrt(epsOxygen * epsChlorine));
        potentialQOCl = new P2Ewald1Real(chargeOxygen * chargeChlorine, ewaldParams.alpha);
        potentialMaster.setPairPotential(oType, clType, tf.make(potentialLJOCl, potentialQOCl));

        potentialLJONa = new P2LennardJones(space, 0.5 * (sigOxygen + sigSodium), Math.sqrt(epsOxygen * epsSodium));
        potentialQONa = new P2Ewald1Real(chargeOxygen * chargeSodium, ewaldParams.alpha);
        potentialMaster.setPairPotential(oType, naType, tf.make(potentialLJONa, potentialQONa));

        potentialQHNa = new P2Ewald1Real(chargeHydrogen * chargeSodium, ewaldParams.alpha);
        potentialMaster.setPairPotential(hType, naType, tf.make(potentialQHNa));

        potentialQHCl = new P2Ewald1Real(chargeHydrogen * chargeChlorine, ewaldParams.alpha);
        potentialMaster.setPairPotential(hType, clType, tf.make(potentialQHCl));


        potentialMM = new P2LennardJones(space, sigMembrane, epsMembrane);
        potentialMaster.setPairPotential(mType, mType, tf.make(potentialMM));

        potentialMO = new P2LennardJones(space, 0.5 * (sigMembrane + sigOxygen), Math.sqrt(epsMembrane * epsOxygen));
        potentialMaster.setPairPotential(mType, oType, tf.make(potentialMO));

        potentialMNa = new P2LennardJones(space, 0.5 * (sigMembrane + sigSodium), Math.sqrt(epsMembrane * epsSodium));
        potentialMaster.setPairPotential(mType, naType, tf.make(potentialMNa));

        potentialMCl = new P2LennardJones(space, 0.5 * (sigMembrane + sigChlorine), Math.sqrt(epsMembrane * epsChlorine));
        potentialMaster.setPairPotential(mType, clType, tf.make(potentialMCl));

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
        PotentialComputeAggregate pcAgg = new PotentialComputeAggregate(potentialMaster, ewaldFourier, pmBonding, pcField);

        //controller and integrator
        integrator = new IntegratorVelocityVerlet(pcAgg, random, 0.002, Kelvin.UNIT.toSim(298), box);
        integrator.setIsothermal(true);
        integrator.setThermostatInterval(100);
        integrator.setIsothermal(false);
        double blOH = ConformationWater3P.bondLengthOH;
        double angleHOH = ConformationWater3P.angleHOH;
        double blHH = blOH * Math.sqrt(2 - 2 * Math.cos(angleHOH));
        SpeciesAgentManager<ShakeListener.BondConstraints> bondConstraintsManager = new SpeciesAgentManager<>(new SpeciesAgentManager.AgentSource<ShakeListener.BondConstraints>() {
            @Override
            public ShakeListener.BondConstraints makeAgent(ISpecies type) {
                if (type == speciesSolvent) {
                    return new ShakeListener.BondConstraints(new int[][]{{0, 1}, {0, 2}, {1, 2}}, new double[]{blHH, blOH, blOH});
                }
                return null;
            }

            @Override
            public void releaseAgent(ShakeListener.BondConstraints agent, ISpecies type) {
            }
        }, getSpeciesManager());
        ShakeListener shake = new ShakeListener(getSpeciesManager(), bondConstraintsManager, integrator);
        shake.setMaxIterations(100);
        shake.setShakeTolerance(1e-9);
        integrator.getEventManager().addListener(shake);
        getController().addActivity(new ActivityIntegrate(integrator));

        potentialTether = new P1Tether(box, speciesMembrane, space);
        potentialTether.setEpsilon(20000 * 298.0 / 125);
        pcField.setFieldPotential(speciesMembrane.getLeafType(), potentialTether);

        BoxImposePbc pbc = new BoxImposePbc(box, space);
        pbc.setApplyToMolecules(true);
        integrator.getEventManager().addListener(new IntegratorListenerAction(pbc));
    }

    public double getWallForce() {
        return wallForce;
    }

    public static void main(String[] args) {
        ReverseOsmosisWater sim = new ReverseOsmosisWater();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, Long.MAX_VALUE));
    }//end of main

}
