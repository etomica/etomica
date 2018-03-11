package etomica.parser.schroeder;

import com.fasterxml.jackson.core.JsonParser;
import com.fasterxml.jackson.databind.ObjectMapper;
import etomica.action.SimulationRestart;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeKinetic;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.*;
import etomica.space2d.Space2D;
import etomica.species.Species;
import etomica.species.SpeciesSpheres;
import etomica.species.SpeciesSpheresMono;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

public class SchroederMDSim extends Simulation {
    private final MDConfig config;
    private final IntegratorVelocityVerlet integrator;
    public final Configuration configuration;

    public SchroederMDSim(MDConfig config) {
        super(Space2D.getInstance());
        this.config = config;

        SpeciesSpheresMono species = new SpeciesSpheresMono(this.space, AtomType.simple("A", 1.0));
        species.setIsDynamic(true);
        SpeciesSpheresMono fixed = new SpeciesSpheresMono(this.space, AtomType.simple("F", Double.POSITIVE_INFINITY));
        fixed.setIsDynamic(true);
        this.addSpecies(species);
        this.addSpecies(fixed);

        Boundary boundary = new BoundaryRectangularNonperiodic(this.space, config.size);
//        Boundary boundary = new BoundaryRectangularPeriodic(this.space, config.size);
        Box box = this.makeBox(boundary);

//        box.setNMolecules(species, config.N);
//        box.setNMolecules(fixed, config.fixedList.length);
        configuration = b -> {
            Set<Integer> fixedIndices = new HashSet<>();
            Arrays.stream(config.fixedList).forEach(fixedIndices::add);
            for (int i = 0; i < config.data.length; i++) {
                double[] atomData = config.data[i];
                IMolecule mol = fixedIndices.contains(i) ? fixed.makeMolecule() : species.makeMolecule();
                IAtomKinetic atom = (IAtomKinetic) mol.getChildList().get(0);
                atom.getPosition().E(Vector.of(atomData[0] - (config.size / 2), -atomData[1] + (config.size / 2)));
                atom.getVelocity().E(Vector.of(atomData[2], atomData[3]));
                b.addMolecule(mol);
            }
        };

        configuration.initializeCoordinates(box);

//        PotentialMasterList pm = new PotentialMasterList(this, this.space);
        PotentialMaster pm = new PotentialMasterMonatomic(this);

        P2SoftSphericalTruncated p2 = new P2SoftSphericalTruncated(this.space, new P2LennardJones(this.space), 3);
        P1WCAWall pWallX = new P1WCAWall(this.space, 0);
        P1WCAWall pWallY = new P1WCAWall(this.space, 1);
        pWallX.setBox(box);
        pWallY.setBox(box);

        pm.addPotential(p2, new AtomType[]{species.getLeafType(), species.getLeafType()});
        pm.addPotential(p2, new AtomType[]{species.getLeafType(), fixed.getLeafType()});
        pm.addPotential(pWallX, new AtomType[]{species.getLeafType()});
        pm.addPotential(pWallY, new AtomType[]{species.getLeafType()});

        integrator = new IntegratorVelocityVerlet(this, pm, box);
        integrator.setTimeStep(config.dt);
//        integrator.getEventManager().addListener(pm.getNeighborManager(box));
    }

    public static void main(String[] args) throws IOException {
        ObjectMapper om = new ObjectMapper();
        om.configure(JsonParser.Feature.ALLOW_UNQUOTED_FIELD_NAMES, true);
        MDConfig config = om.readValue(new File("md.json"), MDConfig.class);

        SchroederMDSim sim = new SchroederMDSim(config);
        ActivityIntegrate ai = new ActivityIntegrate(sim.integrator);
        sim.getController().addAction(ai);
        SimulationGraphic graphic = new SimulationGraphic(sim);
//        SimulationRestart restart = new SimulationRestart(sim);
//        restart.setConfiguration(sim.configuration);
//        graphic.getController().getReinitButton().setAction(restart);
        graphic.getController().getReinitButton().setPostAction(graphic.getPaintAction(sim.box()));
        graphic.makeAndDisplayFrame();
    }
}
