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
import etomica.integrator.IntegratorMD;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.LatticeCubicFcc;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.SpeciesSpheresMono;
import etomica.tests.TestLJMC3D;
import etomica.util.random.RandomMersenneTwister;

import java.foreign.Scope;

public class LJMD3DVecSys extends Simulation {
//    public IntegratorVelocityVerletLessFast integrator;
    public IntegratorMD integrator;
    public SpeciesSpheresMono species;
    public Box box;
    public P2LennardJones potential;
    public Controller controller;
    public DataSourceScalar energy;
    public AccumulatorAverageCollapsing avgEnergy;
    public DataPump pump;

    public LJMD3DVecSys(String type, Scope scope, int nAtoms) {
        super(Space3D.getInstance());
        this.setRandom(new RandomMersenneTwister(1234));
        System.out.println("Using " + type);

        species = new SpeciesSpheresMono(this, space);
        species.setIsDynamic(true);
        addSpecies(species);

        PotentialMaster potentialMaster = new PotentialMasterMonatomic(this);
        box = this.makeBox();
        box.setNMolecules(species, nAtoms);
//        ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(space), space);

        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(0.65);
        inflater.actionPerformed();
//        Configuration configuration = Configurations.fromResourceFile(String.format("LJMC3D%d.pos", 4000), TestLJMC3D.class);
//        configuration.initializeCoordinates(box);
        ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        configuration.initializeCoordinates(box);

        var lj = new P2LennardJones(space, 1.0, 1.0);
        var trunc = new P2SoftSphericalTruncated(lj, 3);
        var leafType = species.getLeafType();
        potentialMaster.addPotential(trunc, new AtomType[]{leafType, leafType});

        switch (type) {
            case "baseline":
                integrator = new IntegratorVelocityVerlet(potentialMaster, this.getRandom(), 0.02, 1, box);
                energy = new MeterPotentialEnergy(potentialMaster, box);
                break;
            case "objects":
                integrator = new IntegratorVelocityVerletLessFast(potentialMaster, this.getRandom(), 0.02, 1, box);
                energy = ((EnergyMeter) integrator).getMeter();
                break;
            case "vecsys1":
                integrator = new IntegratorVelocityVerletFast(potentialMaster, this.getRandom(), 0.02, 1, box);
                energy = ((EnergyMeter) integrator).getMeter();
                break;
            case "vecsys2":
                integrator = new IVVSimd(potentialMaster, this.getRandom(), 0.02, 1, box);
                energy = ((EnergyMeter) integrator).getMeter();
                break;
            case "native":
                integrator = new ComputeForcesForeign(potentialMaster, this.getRandom(), 0.02, 1, box, scope);
                energy = ((ComputeForcesForeign) integrator).getMeter();
                break;
        }

        avgEnergy = new AccumulatorAverageCollapsing();
        avgEnergy.setPushInterval(10);
        pump = new DataPump(energy, avgEnergy);
        IntegratorListenerAction pumpListener = new IntegratorListenerAction(pump);
        pumpListener.setInterval(10);
        integrator.getEventManager().addListener(pumpListener);
    }

    public static void main(String[] args) {
        try (Scope sc = Scope.newNativeScope()) {
            LJMD3DVecSys sim = new LJMD3DVecSys("baseline", sc, 4000);
//            LJMD3DVecSys sim = new LJMD3DVecSys("objects", sc);
//            LJMD3DVecSys sim = new LJMD3DVecSys("vecsys1", sc);
            ActionIntegrate ai = new ActionIntegrate(sim.integrator);
            ai.setMaxSteps(100);
            sim.getController().addAction(ai);

            long t0 = System.nanoTime();
            sim.getController().actionPerformed();
            long t1 = System.nanoTime();

            System.out.println((t1 - t0) / 1_000_000);
        }
    }
}
