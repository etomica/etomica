package etomica.integrationtests;

import com.fasterxml.jackson.databind.ObjectMapper;
import etomica.action.ActionIntegrate;
import etomica.action.BoxInflate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorHard;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2HardSphere;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.random.RandomMersenneTwister;
import org.junit.Test;

import java.io.IOException;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

import static org.junit.Assert.assertEquals;

public class TestSnapshots {
    private static final ObjectMapper om = new ObjectMapper();

    @Test
    public void testHSMD3DNeighborListSnapshot() throws URISyntaxException, IOException {
        String coordsStr = Files.lines(Paths.get(TestSnapshots.class.getResource("HSMD3DNeighborList_021af6881.json").toURI()))
                .findFirst().get();

        Simulation sim = new HSMD3DNeighborList();
        sim.getController().actionPerformed();
        List<Vector> coords = sim.box().getLeafList().getAtoms().stream()
                .map(IAtom::getPosition)
                .collect(Collectors.toList());

        assertEquals(coordsStr, om.writeValueAsString(coords));
    }

    private static class HSMD3DNeighborList extends Simulation {
        public HSMD3DNeighborList() {
            super(Space3D.getInstance());
            this.setRandom(new RandomMersenneTwister(new int[]{1, 2, 3, 4}));

            SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
            species.setIsDynamic(true);
            addSpecies(species);

            addBox(new Box(this.space));

            PotentialMasterList pm = new PotentialMasterList(this, 1.6, this.space);
            P2HardSphere potential = new P2HardSphere(space, 1.0, true);
            AtomType leafType = species.getLeafType();

            pm.addPotential(potential, new AtomType[]{leafType, leafType});

            box().setNMolecules(species, 256);
            BoxInflate inflater = new BoxInflate(box(), space);
            inflater.setTargetDensity(.35 * 2 * space.D() / Math.PI);
            inflater.actionPerformed();
            new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box());

            IntegratorHard integrator = new IntegratorHard(this, pm, box());
            integrator.setIsothermal(false);
            integrator.setTimeStep(0.01);

            ActionIntegrate ai = new ActionIntegrate(integrator);
            ai.setMaxSteps(500);
            getController().addAction(ai);

            integrator.getEventManager().addListener(pm.getNeighborManager(box()));
        }

    }

//    public static void main(String[] args) throws IOException {
//        Simulation sim = new HSMD3DNeighborList();
//
//        List<Vector> coords = sim.box().getLeafList().getAtoms().stream()
//                .map(IAtom::getPosition)
//                .collect(Collectors.toList());
//
//        ObjectMapper om = new ObjectMapper();
//        System.out.println(om.writerWithDefaultPrettyPrinter().writeValueAsString(coords));
//
//        sim.getController().actionPerformed();
//
//        List<Vector> coords2 = sim.box().getLeafList().getAtoms().stream()
//                .map(IAtom::getPosition)
//                .collect(Collectors.toList());
//
//        Process proc = new ProcessBuilder("git", "rev-parse", "--short", "HEAD").start();
//        String hash = new BufferedReader(new InputStreamReader(proc.getInputStream())).readLine();
//        System.out.println(hash);
//        File f = new File(String.format("HSMD3DNeighborList_%s.json", hash));
//        om.writeValue(f, coords2);
//        System.out.println(om.writerWithDefaultPrettyPrinter().writeValueAsString(coords2));
//    }
}
