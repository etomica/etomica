package etomica.modules.joulethomson;

import etomica.action.Action;
import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.graphics.DisplayPhase;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorGear4NPH;
import etomica.integrator.IntegratorMD;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.integrator.IntervalActionAdapter;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.lattice.SpaceLattice;
import etomica.phase.Phase;
import etomica.potential.P2LennardJones;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Bar;
import etomica.units.CompoundUnit;
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
    
    public JouleThomsonSim() {this(Space2D.getInstance());}
    public JouleThomsonSim(Space space) {
        super(space);
        int nAtoms = (space.D() < 3) ? 50 : 64;
        
		defaults.historyPeriod = 1000;
 
        //integrator
        integratorNVE = new IntegratorVelocityVerlet(this);
        integrator = new IntegratorGear4NPH(this);
        integrator.setRelaxationRateP(500.);
        integrator.setRelaxationRateH(300.);
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
	    potential = new P2LennardJones(this);
        potentialMaster.addPotential(potential, new Species[]{species, species});
	    phase = new Phase(this);
        phase.getAgent(species).setNMolecules(nAtoms);
        
        SpaceLattice lattice;
        if (space.D() == 2) {
            lattice = new LatticeOrthorhombicHexagonal();
        }
        else {
            lattice = new LatticeCubicFcc();
        }
        Configuration config = new ConfigurationLattice(lattice);
        config.initializeCoordinates(phase);
        
        integratorJT = new IntegratorJT(potentialMaster, integrator, integratorNVE);
        integratorJT.addListener(new IntervalActionAdapter(new PhaseImposePbc(phase)));
        integrator.setPhase(phase);
        integratorNVE.setPhase(phase);

        integrator.setTimeStep(0.005);
        integratorNVE.setTimeStep(0.005);

        activityIntegrate = new ActivityIntegrate(this, integratorJT);
        getController().addAction(activityIntegrate);
    }
    
    public static void main(String[] args) {
        JouleThomsonSim sim = new JouleThomsonSim(Space3D.getInstance());
        SimulationGraphic simGraphic = new SimulationGraphic(sim);
        simGraphic.makeAndDisplayFrame();
        final DisplayPhase displayPhase = simGraphic.getDisplayPhase(sim.phase);
        sim.activityIntegrate.setDoSleep(true);
        sim.activityIntegrate.setSleepPeriod(10);

        sim.integratorJT.addListener(new IntervalActionAdapter(new Action() {
            public void actionPerformed() {
                displayPhase.repaint();
            }
            public String getLabel() {return "";}
        }));

    }
}