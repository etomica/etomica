/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.glass;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorMD;
import etomica.integrator.IntegratorMD.ThermostatType;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

public class SimGlass extends Simulation {

    public SpeciesSpheresMono speciesA, speciesB;
    public Box box;
    public IntegratorMD integrator;
    public ActivityIntegrate activityIntegrate;
    public final MCMoveSwap swapMove;
    public final IntegratorMC integratorMC;
    public final PotentialChoice potentialChoice;
    public int log2StepS, log2StepE;
    public double minDrFilter;
    public enum PotentialChoice {LJ, WCA, SS, HS}

    public SimGlass(int D, int nA, int nB, double density, double temperature, boolean doSwap, PotentialChoice pc, int log2StepS, int log2StepE, double minDrFilter) {
        super(Space.getInstance(D));
        this.potentialChoice = pc;
        this.log2StepS = log2StepS;
        this.log2StepE = log2StepE;
        this.minDrFilter = minDrFilter;
        //species
        speciesA = new SpeciesSpheresMono(this, space);
        speciesA.setIsDynamic(true);
        addSpecies(speciesA);
        speciesB = new SpeciesSpheresMono(this, space);
        speciesB.setIsDynamic(true);
        addSpecies(speciesB);

        PotentialMasterList potentialMaster = new PotentialMasterList(this, 2.99, space);

        //controller and integrator
        box = this.makeBox();

        integrator = potentialChoice == PotentialChoice.HS ?
                new IntegratorHard(potentialMaster, random, 0.005, temperature, box) :
                new IntegratorVelocityVerlet(potentialMaster, random, 0.005, temperature, box);
        integrator.setIsothermal(true);
        integrator.setThermostat(ThermostatType.ANDERSEN);
        integrator.setThermostatInterval(1);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        integrator.setThermostatNoDrift(true);

        int chs = 85;
        double coreHS = 0.01 * chs;
        P2SquareWell p2AA = null;
        if (potentialChoice == PotentialChoice.LJ) {

            P2LennardJones potentialAA = new P2LennardJones(space);
            P2SoftSphericalTruncated p2TruncatedAA = new P2SoftSphericalTruncatedForceShifted(space, potentialAA, 2.5);
            potentialMaster.addPotential(p2TruncatedAA, new AtomType[]{speciesA.getLeafType(), speciesA.getLeafType()});
            P2LennardJones potentialAB = new P2LennardJones(space, 0.8, 1.5);
            P2SoftSphericalTruncated p2TruncatedAB = new P2SoftSphericalTruncatedForceShifted(space, potentialAB, 2.5);
            potentialMaster.addPotential(p2TruncatedAB, new AtomType[]{speciesA.getLeafType(), speciesB.getLeafType()});
            P2LennardJones potentialBB = new P2LennardJones(space, 0.88, 0.5);
            P2SoftSphericalTruncated p2TruncatedBB = new P2SoftSphericalTruncatedForceShifted(space, potentialBB, 2.5);
            potentialMaster.addPotential(p2TruncatedBB, new AtomType[]{speciesB.getLeafType(), speciesB.getLeafType()});
        } else if (potentialChoice == PotentialChoice.WCA) {
            potentialMaster.setRange(2);
            // https://doi.org/10.1103/PhysRevX.1.021013
            P2WCA potentialAA = new P2WCA(space, 1, 1);
            potentialMaster.addPotential(potentialAA, new AtomType[]{speciesA.getLeafType(), speciesA.getLeafType()});
            P2WCA potentialAB = new P2WCA(space, 0.5 + 0.5 / 1.4, 1);
            potentialMaster.addPotential(potentialAB, new AtomType[]{speciesA.getLeafType(), speciesB.getLeafType()});
            P2WCA potentialBB = new P2WCA(space, 1.0 / 1.4, 1);
            potentialMaster.addPotential(potentialBB, new AtomType[]{speciesB.getLeafType(), speciesB.getLeafType()});
        } else if (potentialChoice == PotentialChoice.SS) {
            // https://doi.org/10.1103/PhysRevLett.81.120 prescribes cut=4.5*(0.5+0.5/1.4)=3.85714
            P2SoftSphere potentialAA = new P2SoftSphere(space, 1, 1, 12);
            P2SoftSphericalTruncated p2TruncatedAA = new P2SoftSphericalTruncatedForceShifted(space, potentialAA, 2.5);
            potentialMaster.addPotential(p2TruncatedAA, new AtomType[]{speciesA.getLeafType(), speciesA.getLeafType()});
            P2SoftSphere potentialAB = new P2SoftSphere(space, 0.5 + 0.5 / 1.4, 1, 12);
            P2SoftSphericalTruncated p2TruncatedAB = new P2SoftSphericalTruncatedForceShifted(space, potentialAB, 2.5);
            potentialMaster.addPotential(p2TruncatedAB, new AtomType[]{speciesA.getLeafType(), speciesB.getLeafType()});
            P2SoftSphere potentialBB = new P2SoftSphere(space, 1.0 / 1.4, 1, 12);
            P2SoftSphericalTruncated p2TruncatedBB = new P2SoftSphericalTruncatedForceShifted(space, potentialBB, 2.5);
            potentialMaster.addPotential(p2TruncatedBB, new AtomType[]{speciesB.getLeafType(), speciesB.getLeafType()});
        } else if (potentialChoice == PotentialChoice.HS) {
            double L = Math.pow((nA + nB) / density, 1.0 / 3.0);
            if (L < 2.01) throw new RuntimeException("too small!");
            double nbrCut = 1.7;
            if (L < nbrCut * 2) nbrCut = L / 2.001;
            potentialMaster.setRange(nbrCut);
            p2AA = new P2SquareWell(space, coreHS, 1 / coreHS, -100, false);
            potentialMaster.addPotential(p2AA, new AtomType[]{speciesA.getLeafType(), speciesA.getLeafType()});
            P2HardSphere potentialAB = new P2HardSphere(space, 0.5 + 0.5 / 1.4, false);
            potentialMaster.addPotential(potentialAB, new AtomType[]{speciesA.getLeafType(), speciesB.getLeafType()});
            P2HardSphere potentialBB = new P2HardSphere(space, 1.0 / 1.4, false);
            potentialMaster.addPotential(potentialBB, new AtomType[]{speciesB.getLeafType(), speciesB.getLeafType()});
            integrator.setAlwaysScaleRandomizedMomenta(true);
        }

        //construct box
        box.setNMolecules(speciesA, nA);
        box.setNMolecules(speciesB, nB);
        BoxInflate boxInflate = new BoxInflate(box, box.getSpace(), density);
        boxInflate.actionPerformed();
        new ConfigurationLattice(space.D() == 2 ? (new LatticeOrthorhombicHexagonal(space)) : (new LatticeCubicFcc(space)), space).initializeCoordinates(box);

        integrator.getEventManager().addListener(potentialMaster.getNeighborManager(box));

        swapMove = new MCMoveSwap(space, random, potentialMaster, speciesA, speciesB);
        integratorMC = new IntegratorMC(potentialMaster, random, integrator.getTemperature(), box);
        integratorMC.getMoveManager().addMCMove(swapMove);
        if (doSwap) {
            integrator.setThermostat(ThermostatType.HYBRID_MC);
            integrator.setIntegratorMC(integratorMC, 10000);
            integrator.setThermostatInterval(1000);
        }

        integrator.reset();
        integrator.doThermostat();
        if (potentialChoice == PotentialChoice.HS) {
            boolean success = false;
            MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster, box);
            integrator.setTimeStep(0.001);
            for (; chs <= 100; chs++) {
                p2AA.setCoreDiameter(chs * 0.01);
                p2AA.setLambda(1 / (chs * 0.01));
                double u = meterPE.getDataAsScalar();
                if (u == Double.POSITIVE_INFINITY) {
                    chs--;
                    p2AA.setCoreDiameter(chs * 0.01);
                    p2AA.setLambda(1 / (chs * 0.01));
                    break;
                }
                if (chs == 100) {
                    success = true;
                    potentialMaster.removePotential(p2AA);
                    P2HardSphere p = new P2HardSphere(space, 1, false);
                    potentialMaster.addPotential(p, new AtomType[]{speciesA.getLeafType(), speciesA.getLeafType()});
                }
            }
            if (!success) {

                while (chs < 100) {
                    integrator.reset();
                    for (int i = 0; i < 1000; i++) {
                        integrator.doStep();
                    }
                    chs++;
                    p2AA.setCoreDiameter(chs * 0.01);
                    p2AA.setLambda(1 / (chs * 0.01));
                    if (meterPE.getDataAsScalar() == Double.POSITIVE_INFINITY) {
                        chs--;
                        p2AA.setCoreDiameter(chs * 0.01);
                        p2AA.setLambda(1 / (chs * 0.01));
                        continue;
                    } else if (chs == 100) {
                        potentialMaster.removePotential(p2AA);
                        P2HardSphere p = new P2HardSphere(space, 1, false);
                        potentialMaster.addPotential(p, new AtomType[]{speciesA.getLeafType(), speciesA.getLeafType()});
                    }
                }

            }
            integrator.reset();
            integrator.resetStepCount();
            integrator.setTimeStep(0.005);
        }

    }

    public static void main(String[] args) {

        GlassParams params = new GlassParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
        }
        SimGlass sim = new SimGlass(params.D, params.nA, params.nB, params.density, params.temperature, params.doSwap, params.potential, params.log2StepS, params.log2StepE, params.minDrFilter);
        sim.getController().actionPerformed();
    }//end of main

    public static class GlassParams extends ParameterBase {
        public int D = 2;
        public int nA = 130, nB = 70;
        // rho=1000/29.34^2=1.16 (yields P=0 at T=0) for LJ http://dx.doi.org/10.1088/0953-8984/21/3/035117
        //     for 3D, 1000/8.88^3 = 1.43
        //       or (with xA=0.8), rho=1.25 yields T0=1.06 (http://dx.doi.org/10.1103/PhysRevX.1.021013)
        // for SS, P=18.37 https://doi.org/10.1103/PhysRevLett.81.120
        //    rho=1.35 at T=Tg=0.55
        // for HS, rho=1.35 is very high, but 1.30 is perhaps not high enough
        // for WCA, rho=0.75*1.4*1.4 = 1.47
        public double density = 1000 / (29.34 * 29.34);
        public double temperature = 1.0;
        public boolean doSwap = true;
        public PotentialChoice potential = PotentialChoice.LJ;
        public int log2StepS =  1;
        public int log2StepE = 40;
        public double minDrFilter = 0.5;
    }
}