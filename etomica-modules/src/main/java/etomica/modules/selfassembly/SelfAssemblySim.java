/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.selfassembly;

import etomica.action.activity.ActivityIntegrate;
import etomica.action.controller.Controller;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.iterator.ApiBuilder;
import etomica.box.Box;
import etomica.config.ConformationLinear;
import etomica.data.meter.MeterTemperature;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMD.ThermostatType;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.config.ConfigurationLatticeRandom;
import etomica.nbr.CriterionAll;
import etomica.nbr.CriterionBondedSimple;
import etomica.nbr.CriterionInterMolecular;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.species.SpeciesSpheresHetero;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Kelvin;

public class SelfAssemblySim extends Simulation {

    public final PotentialMasterList potentialMaster;
    public final ConfigurationLatticeRandom config;
	public IntegratorHard integratorHard;
	public java.awt.Component display;
	public Box box;
	public MeterTemperature thermometer;
	public SpeciesSpheresMono speciesA;
	public SpeciesSpheresHetero speciesB;
	public final P2SquareWell p2AA, p2AB1, p2AB2, p2B1B1, p2B1B2, p2B2B2;
	public final P2Tether p2Bond;
    public ActivityIntegrate activityIntegrate;
    public AtomLeafAgentManager<IAtom[]> agentManager = null;
    public int nA, nB;
    public int nB1, nB2;
    double sigAA, sigAB1, sigAB2, sigB1B1, sigB1B2, sigB2B2;
    double epsAA, epsAB1, epsAB2, epsB1B1, epsB1B2, epsB2B2;
    double lamAA, lamAB1, lamAB2, lamB1B1, lamB1B2, lamB2B2;
    final public AtomType typeB1 = AtomType.simple("B1",1.0);
    final public AtomType typeB2 = AtomType.simple("B2", 1.0);
    final public AtomType typeA = AtomType.simple("A", 1.0);

    public SelfAssemblySim(Space space) {
        super(space);

        nA = 100;
        nB = 10;
        nB1 = 1;
        nB2 = 5;

        speciesA = new SpeciesSpheresMono(space, typeA);
        speciesA.setIsDynamic(true);
        speciesB = new SpeciesSpheresHetero(space, new AtomType[] {typeB1, typeB2});
        speciesB.setConformation(new ConformationLinear(space,0.501));
        speciesB.setChildCount(new int[] {nB1, nB2});

        speciesB.setIsDynamic(true);
        addSpecies(speciesA);
        addSpecies(speciesB);

        potentialMaster = new PotentialMasterList(this, 3, space);
        potentialMaster.setCellRange(1);

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

        integratorHard = new IntegratorHard(this, potentialMaster, box);
        integratorHard.setIsothermal(true);
        integratorHard.setTemperature(Kelvin.UNIT.toSim(300));
        integratorHard.setTimeStep(0.002);
        integratorHard.setThermostat(ThermostatType.ANDERSEN_SINGLE);
        integratorHard.setThermostatInterval(100);
        integratorHard.getEventManager().addListener(((PotentialMasterList) potentialMaster).getNeighborManager(box));

        //potentials
        p2AA = new P2SquareWell(space, sigAA, lamAA, epsAA, true);
        p2AB1 = new P2SquareWell(space, sigAB1, lamAB1, epsAB1, true);
        p2AB2 = new P2SquareWell(space, sigAB2, lamAB2, epsAB2, true);
        p2B1B1 = new P2SquareWell(space, sigB1B1, lamB1B1, epsB1B1, true);
        p2B1B2 = new P2SquareWell(space, sigB1B2, lamB1B2, epsB1B2, true);
        p2B2B2 = new P2SquareWell(space, sigB2B2, lamB2B2, epsB2B2, true);
        p2Bond = new P2Tether(space, sigAA * 0.8, true);

       // p2Bond = new P2HardBond(space, 0.8, 0.2, true);
        PotentialGroup p1Surfactant = potentialMaster.makePotentialGroup(1);
        p1Surfactant.addPotential(p2Bond, ApiBuilder.makeAdjacentPairIterator());
        potentialMaster.addPotential(p1Surfactant, new ISpecies[]{speciesB});

        potentialMaster.addPotential(p2AA, new AtomType[]{typeA, typeA});
        potentialMaster.addPotential(p2AB1, new AtomType[]{typeA, typeB1});
        potentialMaster.addPotential(p2AB2, new AtomType[]{typeA, typeB2});
        potentialMaster.addPotential(p2B1B2, new AtomType[]{typeB1, typeB2});
        potentialMaster.addPotential(p2B1B1, new AtomType[]{typeB1, typeB1});
        potentialMaster.addPotential(p2B2B2, new AtomType[]{typeB2, typeB2});

        CriterionInterMolecular sqwCriterion = (CriterionInterMolecular) potentialMaster.getCriterion(typeB1, typeB1);
        CriterionBondedSimple nonBondedCriterion = new CriterionBondedSimple(new CriterionAll());
        nonBondedCriterion.setBonded(false);
        sqwCriterion.setIntraMolecularCriterion(nonBondedCriterion);

        sqwCriterion = (CriterionInterMolecular) potentialMaster.getCriterion(typeB1, typeB2);
        nonBondedCriterion = new CriterionBondedSimple(new CriterionAll());
        nonBondedCriterion.setBonded(false);
        sqwCriterion.setIntraMolecularCriterion(nonBondedCriterion);

        sqwCriterion = (CriterionInterMolecular) potentialMaster.getCriterion(typeB2, typeB2);
        nonBondedCriterion = new CriterionBondedSimple(new CriterionAll());
        nonBondedCriterion.setBonded(false);
        sqwCriterion.setIntraMolecularCriterion(nonBondedCriterion);


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
        speciesB.setChildCount(new int[] {nB1, nB2});
        box.setNMolecules(speciesB, 0);
        box.setNMolecules(speciesB, nB);

    }
    public int getNB1() {
        return nB1;
    }

    public void setNB2(int nB2) {
        this.nB2 = nB2;
        speciesB.setChildCount(new int[] {nB1, nB2});
        box.setNMolecules(speciesB, 0);
        box.setNMolecules(speciesB, nB);
    }
    public int getNB2() {
        return nB2;
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
}
