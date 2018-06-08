package etomica.experimental;

import etomica.action.ActionIntegrate;
import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.config.Configurations;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.DataPump;
import etomica.data.DataSourceScalar;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.LatticeCubicFcc;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.tests.TestLJMC3D;

public class LJMD3DVecSys extends Simulation {
    public IntegratorVelocityVerletLessFast integrator;
    public SpeciesSpheresMono species;
    public Box box;
    public P2LennardJones potential;
    public Controller controller;
    public DataSourceScalar energy;
    public AccumulatorAverageCollapsing avgEnergy;
    public DataPump pump;


    public LJMD3DVecSys() {
        super(Space3D.getInstance());

        species = new SpeciesSpheresMono(this, space);
        species.setIsDynamic(true);
        addSpecies(species);

        PotentialMaster potentialMaster = new PotentialMasterMonatomic(this);
        box = this.makeBox();
        box.setNMolecules(species, 500);
//        ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(space), space);

        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(0.65);
        inflater.actionPerformed();
        Configuration configuration = Configurations.fromResourceFile(String.format("LJMC3D%d.pos", 500), TestLJMC3D.class);
        configuration.initializeCoordinates(box);

//        integrator = new IntegratorVelocityVerletFast(potentialMaster, this.getRandom(), 0.02, 1, box);
        integrator = new IntegratorVelocityVerletLessFast(potentialMaster, this.getRandom(), 0.02, 1, box);

        energy = integrator.getMeter();
        avgEnergy = new AccumulatorAverageCollapsing();
        avgEnergy.setPushInterval(10);
        pump = new DataPump(energy, avgEnergy);
        IntegratorListenerAction pumpListener = new IntegratorListenerAction(pump);
        pumpListener.setInterval(10);
        integrator.getEventManager().addListener(pumpListener);
    }

    public static void main(String[] args) {
        LJMD3DVecSys sim = new LJMD3DVecSys();
        ActionIntegrate ai = new ActionIntegrate(sim.integrator);
        ai.setMaxSteps(1000000 / 500);
        sim.getController().addAction(ai);

        long t0 = System.nanoTime();
        sim.getController().actionPerformed();
        long t1 = System.nanoTime();

        System.out.println((t1 - t0) / 1_000_000);
    }
}
