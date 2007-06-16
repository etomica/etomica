package etomica.modules.joulethomson;

import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomTypeSphere;
import etomica.chem.elements.ElementSimple;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.graphics.DisplayPhase;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorGear4NPH;
import etomica.integrator.IntegratorMD;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.lattice.SpaceLattice;
import etomica.phase.Phase;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Bar;
import etomica.units.CompoundUnit;
import etomica.units.Kelvin;
import etomica.units.Meter;
import etomica.units.Prefix;
import etomica.units.PrefixedUnit;
import etomica.units.Unit;

public class JouleThomsonSim extends Simulation {
    
    private static final long serialVersionUID = 1L;
    IntegratorGear4NPH integrator;
    IntegratorMD integratorNVE;
    SpeciesSpheresMono species;
    P2LennardJones potential;
    Phase phase;
    IntegratorJT integratorJT;
    ActivityIntegrate activityIntegrate;
    Configuration config;

    public JouleThomsonSim() {this(Space2D.getInstance());}
    public JouleThomsonSim(Space space) {
        super(space);
        PotentialMaster potentialMaster = new PotentialMaster(space);
        int nAtoms = (space.D() < 3) ? 50 : 64;
        double sigma = 3.0;
        
        //integrator
        integratorNVE = new IntegratorVelocityVerlet(this, potentialMaster);
        integrator = new IntegratorGear4NPH(this, potentialMaster);
        integrator.setRelaxationRateP(500.);
        integrator.setRelaxationRateH(300.);
        integratorNVE.setTemperature(Kelvin.UNIT.toSim(300));
        integrator.setTemperature(Kelvin.UNIT.toSim(300));
        final Unit pUnit;
        if (space.D() == 2) {
            Unit[] units = new Unit[] {Bar.UNIT, new PrefixedUnit(Prefix.NANO, Meter.UNIT)};
            double[] exponents = new double[] {1.0, 1.0};
            pUnit = new CompoundUnit(units, exponents);
        }
        else {
            pUnit = Bar.UNIT;
        }
        integrator.setTargetP(pUnit.toSim(100.0));
	    
	    //species and potential
	    species = new SpeciesSpheresMono(this);
        ((AtomTypeSphere)species.getMoleculeType()).setDiameter(sigma);
        ((ElementSimple)((AtomTypeSphere)species.getMoleculeType()).getElement()).setMass(40);
        getSpeciesManager().addSpecies(species);
	    potential = new P2LennardJones(space, sigma, Kelvin.UNIT.toSim(300));
        potentialMaster.addPotential(potential, new Species[]{species, species});
	    phase = new Phase(this);
        addPhase(phase);
        phase.getAgent(species).setNMolecules(nAtoms);
        
        SpaceLattice lattice;
        if (space.D() == 2) {
            lattice = new LatticeOrthorhombicHexagonal();
        }
        else {
            lattice = new LatticeCubicFcc();
        }
        config = new ConfigurationLattice(lattice);
        config.initializeCoordinates(phase);
        
        integratorJT = new IntegratorJT(getRandom(), integrator, integratorNVE);
        integratorJT.addIntervalAction(new PhaseImposePbc(phase));
        integrator.setPhase(phase);
        integratorNVE.setPhase(phase);

        integrator.setTimeStep(0.001);
        integratorNVE.setTimeStep(0.005);

        activityIntegrate = new ActivityIntegrate(integratorJT);
        getController().addAction(activityIntegrate);
    }
    
    public static void main(String[] args) {
        JouleThomsonSim sim = new JouleThomsonSim(Space3D.getInstance());
        SimulationGraphic simGraphic = new SimulationGraphic(sim);
        simGraphic.makeAndDisplayFrame();
        final DisplayPhase displayPhase = simGraphic.getDisplayPhase(sim.phase);
        sim.activityIntegrate.setDoSleep(true);
        sim.activityIntegrate.setSleepPeriod(10);

        sim.integratorJT.addIntervalAction(simGraphic.getDisplayPhasePaintAction(sim.phase));

    }
}