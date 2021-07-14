/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.selfassembly;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.config.ConfigurationLatticeRandom;
import etomica.config.ConformationLinear;
import etomica.data.meter.MeterTemperature;
import etomica.integrator.IntegratorHardFasterer;
import etomica.integrator.IntegratorMDFasterer;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.nbr.list.NeighborListManagerFastererHard;
import etomica.potential.P2HardGeneric;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.compute.PotentialComputePairGeneral;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;
import etomica.species.SpeciesGeneralMutable;
import etomica.units.Kelvin;

import java.util.ArrayList;
import java.util.List;

public class SelfAssemblySimFasterer extends Simulation {

    public final PotentialMasterBonding pcBonding;
    public final ConfigurationLatticeRandom config;
    public IntegratorHardFasterer integratorHard;
    public java.awt.Component display;
    public Box box;
    public MeterTemperature thermometer;
    public SpeciesGeneral speciesA;
    public SpeciesGeneralMutable speciesB;
    public NeighborListManagerFastererHard neighborManager;
    public final P2HardGeneric p2AA, p2AB1, p2AB2, p2B1B1, p2B1B2, p2B2B2;
    public final P2HardGeneric p2Bond;
    public ActivityIntegrate activityIntegrate;
    public AtomLeafAgentManager<IAtom[]> agentManager = null;
    public int nA, nB;
    public int nB1, nB2;
    double sigAA, sigAB1, sigAB2, sigB1B1, sigB1B2, sigB2B2;
    double epsAA, epsAB1, epsAB2, epsB1B1, epsB1B2, epsB2B2;
    double lamAA, lamAB1, lamAB2, lamB1B1, lamB1B2, lamB2B2;
    final public AtomType typeB1 = AtomType.simple("B1", 1.0);
    final public AtomType typeB2 = AtomType.simple("B2", 1.0);
    final public AtomType typeA = AtomType.simple("A", 1.0);

    public SelfAssemblySimFasterer(Space space) {
        super(space);

        nA = 100;
        nB = 10;
        nB1 = 1;
        nB2 = 5;

        speciesA = SpeciesGeneral.monatomic(space, typeA, true);

        speciesB = new SpeciesGeneralMutable(this.makeSpeciesB(nB1, nB2));
        addSpecies(speciesA);
        addSpecies(speciesB);

        double defaultEps = 50;
        sigAA = 1.0;
        epsAA = Kelvin.UNIT.toSim(defaultEps);

        sigAB1 = 1.0;
        epsAB1 = Kelvin.UNIT.toSim(defaultEps);

        sigAB2 = 1.0;
        epsAB2 = Kelvin.UNIT.toSim(defaultEps);

        sigB1B1 = 1.5;
        epsB1B1 = Kelvin.UNIT.toSim(defaultEps);

        sigB1B2 = 1.0;
        epsB1B2 = Kelvin.UNIT.toSim(0.);

        sigB2B2 = 1.0;
        epsB2B2 = Kelvin.UNIT.toSim(defaultEps);

        lamAA = lamAB1 = lamAB2 = lamB1B1 = lamB1B2 = lamB2B2 = 1.5;

        box = this.makeBox(new BoundaryRectangularPeriodic(space, space.D() == 2 ? 60 : 20));
        box.setNMolecules(speciesA, nA);
        box.setNMolecules(speciesB, nB);
        config = new ConfigurationLatticeRandom(space.D() == 2 ? new LatticeOrthorhombicHexagonal(space) : new LatticeCubicFcc(space), space, random);
        config.initializeCoordinates(box);

        pcBonding = new PotentialMasterBonding(getSpeciesManager(), box);

        neighborManager = new NeighborListManagerFastererHard(getSpeciesManager(), box, 2, 3, pcBonding.getBondingInfo());
        PotentialComputePairGeneral pcPair = new PotentialComputePairGeneral(getSpeciesManager(), box, neighborManager);
        //potentials
        p2AA = new P2HardGeneric(new double[]{sigAA, sigAA * lamAA}, new double[]{Double.POSITIVE_INFINITY, -epsAA}, true);
        p2AB1 = new P2HardGeneric(new double[]{sigAB1, sigAB1 * lamAB1}, new double[]{Double.POSITIVE_INFINITY, -epsAB1}, true);
        p2AB2 = new P2HardGeneric(new double[]{sigAB2, sigAB2 * lamAB2}, new double[]{Double.POSITIVE_INFINITY, -epsAB2}, true);
        p2B1B1 = new P2HardGeneric(new double[]{sigB1B1, sigB1B1 * lamB1B1}, new double[]{Double.POSITIVE_INFINITY, -epsB1B1}, true);
        p2B1B2 = new P2HardGeneric(new double[]{sigB1B2, sigB1B2 * lamB1B2}, new double[]{Double.POSITIVE_INFINITY, -epsB1B2}, true);
        p2B2B2 = new P2HardGeneric(new double[]{sigB2B2, sigB2B2 * lamB2B2}, new double[]{Double.POSITIVE_INFINITY, -epsB2B2}, true);
        p2Bond = new P2HardGeneric(new double[]{sigAA * 0.8, sigAA}, new double[]{0, Double.POSITIVE_INFINITY}, false);

        setupBonding();

        pcPair.setPairPotential(typeA, typeA, p2AA);
        pcPair.setPairPotential(typeA, typeB1, p2AB1);
        pcPair.setPairPotential(typeA, typeB2, p2AB2);
        pcPair.setPairPotential(typeB1, typeB2, p2B1B2);
        pcPair.setPairPotential(typeB1, typeB1, p2B1B1);
        pcPair.setPairPotential(typeB2, typeB2, p2B2B2);

        integratorHard = new IntegratorHardFasterer(IntegratorHardFasterer.extractHardPotentials(pcPair), neighborManager, random, 0.002, Kelvin.UNIT.toSim(300), box, getSpeciesManager(), pcBonding.getBondingInfo());
        integratorHard.setIsothermal(true);
        integratorHard.setThermostat(IntegratorMDFasterer.ThermostatType.ANDERSEN_SINGLE);
        integratorHard.setThermostatInterval(100);

        // **** Setting Up the thermometer Meter *****

        thermometer = new MeterTemperature(box, space.D());

        activityIntegrate = new ActivityIntegrate(integratorHard, Long.MAX_VALUE, true);
        getController().addActivity(activityIntegrate);
    }

    public int getNA() {
        return nA;
    }

    public void setNA(int nA) {
        this.nA = nA;
        box.setNMolecules(speciesA, nA);
    }

    public int getNB() {
        return nB;
    }

    public void setNB(int nB) {
        this.nB = nB;
        box.setNMolecules(speciesB, nB);
    }

    public void setNB1(int nB1) {
        this.nB1 = nB1;
        speciesB.setSpecies(this.makeSpeciesB(nB1, nB2));
        box.setNMolecules(speciesB, 0);
        box.setNMolecules(speciesB, nB);
        setupBonding();
    }

    public int getNB1() {
        return nB1;
    }

    public void setNB2(int nB2) {
        this.nB2 = nB2;
        speciesB.setSpecies(this.makeSpeciesB(nB1, nB2));
        box.setNMolecules(speciesB, 0);
        box.setNMolecules(speciesB, nB);
        setupBonding();
    }

    public int getNB2() {
        return nB2;
    }

    public void setupBonding() {
        List<int[]> bondList = new ArrayList<>();
        for (int i = 0; i < nB1 + nB2 - 1; i++) {
            bondList.add(new int[]{i, i + 1});
        }
        pcBonding.setBondingPotentialPair(speciesB, p2Bond, bondList);
    }

    public AtomType getTypeA() {
        return typeA;
    }

    public AtomType getTypeB1() {
        return typeB1;
    }

    public AtomType getTypeB2() {
        return typeB2;
    }

    public SpeciesGeneral makeSpeciesB(int nB1, int nB2) {
        return new SpeciesBuilder(space)
                .addCount(typeB1, nB1)
                .addCount(typeB2, nB2)
                .withConformation(new ConformationLinear(space, 0.501))
                .setDynamic(true)
                .build();
    }
}
