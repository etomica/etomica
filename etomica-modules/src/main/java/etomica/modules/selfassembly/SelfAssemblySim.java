/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.selfassembly;

import etomica.action.activity.ActivityIntegrate;
import etomica.action.controller.Controller;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.meter.MeterTemperature;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMD.ThermostatType;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.molecule.IMolecule;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2HardSphere;
import etomica.potential.P2SquareWell;
import etomica.potential.P2Tether;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresHetero;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Kelvin;

public class SelfAssemblySim extends Simulation {

    public final PotentialMaster potentialMaster;
    public final ConfigurationLatticeRandom config;
    public Controller controller1;
	public IntegratorHard integratorHard;
	public java.awt.Component display;
	public Box box;
	public MeterTemperature thermometer;
	public SpeciesSpheresMono speciesA;
	public SpeciesSpheresHetero speciesB;
	public final P2SquareWell p2AA, p2AB1, p2AB2, p2B1B1, p2B1B2, p2B2B2;
	public final P2Tether ABbonded;
    public ActivityIntegrate activityIntegrate;
    public AtomLeafAgentManager<IAtom[]> agentManager = null;
    public int nA, nB;
    public int nB1, nB2;
    double sigAA, sigAB1, sigAB2, sigB1B1, sigB1B2, sigB2B2;
    double epsAA, epsAB1, epsAB2, epsB1B1, epsB1B2, epsB2B2;
    double lamAA, lamAB1, lamAB2, lamB1B1, lamB1B2, lamB2B2;
    final AtomType typeB1 = AtomType.simple("B1",1.0);
    final AtomType typeB2 = AtomType.simple("B2", 1.0);
    final AtomType typeA = AtomType.simple("A", 1.0);

    public SelfAssemblySim(Space space) {
        super(space);

        speciesA = new SpeciesSpheresMono(space, typeA);
        speciesA.setIsDynamic(true);
        speciesB = new SpeciesSpheresHetero(space, new AtomType[] {typeB1, typeB2});
        speciesB.setIsDynamic(true);
        addSpecies(speciesA);
        addSpecies(speciesB);

        potentialMaster = new PotentialMasterList(this, 3, space);
        ((PotentialMasterList) potentialMaster).setCellRange(1);

        controller1 = getController();

        sigAA = sigAB1 = sigAB2 = sigB1B1 = sigB1B2 = sigB2B2 = 1.0;
        epsAA = epsAB1 = epsAB2 = epsB1B1 =  epsB1B2 = epsB2B2 = 1.0;
        lamAA = lamAB1 = lamAB2 = lamB1B1 =  lamB1B2 = lamB2B2 = 2.0;

        box = this.makeBox(new BoundaryRectangularPeriodic(space, space.D() == 2 ? 60 : 20));
        box.setNMolecules(speciesA, 50);
        nA = 100;
        box.setNMolecules(speciesB, 100);
        nB = 10;
        nB1 = 20;
        nB2 = 20;
        config = new ConfigurationLatticeRandom(space.D() == 2 ? new LatticeOrthorhombicHexagonal(space) : new LatticeCubicFcc(space), space, random);
        config.initializeCoordinates(box);

        integratorHard = new IntegratorHard(this, potentialMaster, box);
        integratorHard.setIsothermal(true);
        integratorHard.setTemperature(Kelvin.UNIT.toSim(300));
        integratorHard.setTimeStep(0.002);
        integratorHard.setThermostat(ThermostatType.ANDERSEN_SINGLE);
        integratorHard.setThermostatInterval(1);
        integratorHard.getEventManager().addListener(((PotentialMasterList) potentialMaster).getNeighborManager(box));

        //potentials
        p2AA = new P2SquareWell(space, sigAA, epsAA, lamAA, true);
        p2AB1 = new P2SquareWell(space, sigAB1, epsAB1, lamAB1, true);
        p2AB2 = new P2SquareWell(space, sigAB2, epsAB2, lamAB2, true);
        p2B1B1 = new P2SquareWell(space, sigB1B1, epsB1B1, lamB1B1, true);
        p2B1B2 = new P2SquareWell(space, sigB1B2, epsB1B2, lamB1B2, true);
        p2B2B2 = new P2SquareWell(space, sigB2B2, epsB2B2, lamB2B2, true);
        ABbonded = new P2Tether(space, sigAA * 0.8, true);

        potentialMaster.addPotential(p2AA, new AtomType[]{typeA, typeA});
        potentialMaster.addPotential(p2AB1, new AtomType[]{typeA, typeB1});
        potentialMaster.addPotential(p2AB2, new AtomType[]{typeA, typeB2});
        potentialMaster.addPotential(p2B1B2, new AtomType[]{typeB1, typeB2});
        potentialMaster.addPotential(p2B1B1, new AtomType[]{typeB1, typeB1});
        potentialMaster.addPotential(p2B2B2, new AtomType[]{typeB2, typeB2});

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
        //do something
    }
    public int getNB1() {
        return nB1;
    }

    public void setNB2(int nB2) {
        this.nB2 = nB2;
        //do something
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
