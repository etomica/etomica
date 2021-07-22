package etomica.integrationtests;

import com.fasterxml.jackson.databind.ObjectMapper;
import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
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
import etomica.species.SpeciesGeneral;
import etomica.util.random.RandomMersenneTwister;
import org.junit.jupiter.api.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class TestSnapshots {
    private static final ObjectMapper om = new ObjectMapper();

    @Test
    public void testHSMD3DNeighborListSnapshot() throws URISyntaxException, IOException {
        String coordsStr = Files.lines(Paths.get(TestSnapshots.class.getResource("HSMD3DNeighborList_021af6881.json").toURI()))
                .findFirst().get();

        HSMD3DNeighborList sim = new HSMD3DNeighborList();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, 500));
        List<Vector> coords = sim.box().getLeafList().getAtoms().stream()
                .map(IAtom::getPosition)
                .collect(Collectors.toList());

        assertEquals(coordsStr, om.writeValueAsString(coords));
    }

    private static class HSMD3DNeighborList extends Simulation {
        public final IntegratorHard integrator;
        public HSMD3DNeighborList() {
            super(Space3D.getInstance());
            this.setRandom(new RandomMersenneTwister(new int[]{1, 2, 3, 4}));

            SpeciesGeneral species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
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

            integrator = new IntegratorHard(random, pm, box());
            integrator.setIsothermal(false);
            integrator.setTimeStep(0.01);

            integrator.getEventManager().addListener(pm.getNeighborManager(box()));
        }

    }

    public static void main(String[] args) throws IOException {
        if (true) return;
        HSMD3DNeighborList sim = new HSMD3DNeighborList();

        List<Vector> coords = sim.box().getLeafList().getAtoms().stream()
                .map(IAtom::getPosition)
                .collect(Collectors.toList());

        ObjectMapper om = new ObjectMapper();
        System.out.println(om.writerWithDefaultPrettyPrinter().writeValueAsString(coords));

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, 500));

        List<Vector> coords2 = sim.box().getLeafList().getAtoms().stream()
                .map(IAtom::getPosition)
                .collect(Collectors.toList());

        Process proc = new ProcessBuilder("git", "rev-parse", "--short", "HEAD").start();
        String hash = new BufferedReader(new InputStreamReader(proc.getInputStream())).readLine();
        System.out.println(hash);
        File f = new File(String.format("HSMD3DNeighborList_%s.json", hash));
        om.writeValue(f, coords2);
        System.out.println(om.writerWithDefaultPrettyPrinter().writeValueAsString(coords2));
    }
}
