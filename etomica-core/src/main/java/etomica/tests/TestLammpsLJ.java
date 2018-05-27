package etomica.tests;

import etomica.action.ActionIntegrate;
import etomica.action.BoxInflate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;

/**
 * Simulation that attempts to replicate a LAMMPS benchmark. The input file is
 * <pre>
 * # 3d Lennard-Jones melt
 *
 * variable	x index 1
 * variable	y index 1
 * variable	z index 1
 *
 * variable	xx equal 20*$x
 * variable	yy equal 20*$y
 * variable	zz equal 20*$z
 *
 * units		lj
 * atom_style	atomic
 *
 * lattice		fcc 0.8442
 * region		box block 0 ${xx} 0 ${yy} 0 ${zz}
 * create_box	1 box
 * create_atoms	1 box
 * mass		1 1.0
 *
 * velocity	all create 1.44 87287 loop geom
 *
 * pair_style	lj/cut 2.5
 * pair_coeff	1 1 1.0 1.0 2.5
 *
 * neighbor	0.3 bin
 * neigh_modify	delay 0 every 20 check no
 *
 * fix		1 all nve
 *
 * run		100
 * </pre>
 */
public class TestLammpsLJ extends Simulation {
    public IntegratorVelocityVerlet integrator;

    public TestLammpsLJ() {
        super(Space3D.getInstance());

        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        species.setIsDynamic(true);
        addSpecies(species);

        Box box = this.makeBox();

        P2SoftSphericalTruncated p2 = new P2SoftSphericalTruncated(
                space,
                new P2LennardJones(space, 1.0, 1.0),
                2.5
        );

        box.setNMolecules(species, 32000);
        BoxInflate inflate = new BoxInflate(box, space);
        inflate.setTargetDensity(0.8442);
        inflate.actionPerformed();
        Configuration config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        config.initializeCoordinates(box);

        PotentialMasterList pm = new PotentialMasterList(this, 2.8, space);
        pm.addPotential(p2, new AtomType[]{species.getLeafType(), species.getLeafType()});

//        pm.getNeighborManager(box).setUpdateInterval(2);

        integrator = new IntegratorVelocityVerlet(this, pm, box);
        integrator.setIsothermal(false);
        integrator.setTemperature(1.44);
        integrator.setTimeStep(0.005);
        integrator.getEventManager().addListener(pm.getNeighborManager(box));

    }

    public static void main(String[] args) {
        TestLammpsLJ sim = new TestLammpsLJ();
        ActionIntegrate ai = new ActionIntegrate(sim.integrator);
        ai.setMaxSteps(100);
        sim.getController().addAction(ai);

        long t1 = System.nanoTime();
        sim.getController().actionPerformed();
        long t2 = System.nanoTime();

        System.out.println((t2 - t1) / 1_000_000 + " ms");
    }
}
