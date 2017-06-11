package etomica.freeenergy.npath;

import etomica.box.Box;
import etomica.graphics.DisplayBox;
import etomica.graphics.SimulationGraphic;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

/**
 * This reads ATOM files from LAMMPS and creates a display of one of the
 * configurations.
 * <p>
 * Created by andrew on 4/25/17.
 */
public class ConfigFromFileLAMMPS {

    public static void main(String[] args) throws FileNotFoundException, IOException {
        ConfigFromFileLAMMPSParam params = new ConfigFromFileLAMMPSParam();
        if (args.length == 0) {
            throw new RuntimeException("Usage: ConfigFromFileLAMMPS -filename sim.atom -configNum n");
        }
        ParseArgs.doParseArgs(params, args);
        String filename = params.filename;
        int which = params.configNum;
        FileReader fileReader = new FileReader(filename);
        BufferedReader reader = new BufferedReader(fileReader);
        String line = null;
        String read = "";
        int numAtoms = -1;
        int item = -1;
        Space space = Space3D.getInstance();
        Vector[] edges = space.makeVectorArray(3);
        double[][] boundaryStuff = new double[3][3];
        Simulation sim = new Simulation(space);
        Species species = new SpeciesSpheresMono(sim, space);
        sim.addSpecies(species);
        Box box = new Box(space);
        sim.addBox(box);
        int config = 1;
        while ((line = reader.readLine()) != null) {
            if (line.matches("ITEM:.*")) {
                read = "";
                if (line.equals("ITEM: NUMBER OF ATOMS")) {
                    read = "numAtoms";
                }
                else if (line.matches("ITEM: BOX BOUNDS.*")) {
                    read = "boundary";
                }
                else if (line.matches("ITEM: ATOMS.*")) {
                    read = "atoms";
                    Boundary boundary = new BoundaryDeformablePeriodic(space, edges);
                    box.setBoundary(boundary);
                }
                item = 0;
                continue;
            }
            if (read == "numAtoms") {
                numAtoms = Integer.parseInt(line);
                System.out.println("numAtoms: " + numAtoms);
                box.setNMolecules(species, numAtoms);
            }
            else if (read == "boundary") {
                String[] bits = line.split("[ \t]+");
                if (bits.length == 2) {
                    double xlo = Double.parseDouble(bits[0]);
                    double xhi = Double.parseDouble(bits[1]);
                    edges[item].setX(item, xhi - xlo);
                }
                else {
                    for (int i = 0; i < 3; i++) {
                        boundaryStuff[item][i] = Double.parseDouble(bits[i]);
                    }
                    if (item == 2) {
                        double xy = boundaryStuff[0][2];
                        double xz = boundaryStuff[1][2];
                        double yz = boundaryStuff[2][2];
                        double zlo = boundaryStuff[2][0];
                        double zhi = boundaryStuff[2][1];
                        double yhi = boundaryStuff[1][1] - Math.max(0, yz);
                        double ylo = boundaryStuff[1][0] - Math.min(0, yz);
                        double offset = Math.max(Math.max(Math.max(0, xy), xz), xy + xz);
                        double xhi = boundaryStuff[0][1] - offset;
                        offset = Math.min(Math.min(Math.min(0, xy), xz), xy + xz);
                        double xlo = boundaryStuff[0][0] - offset;
                        edges[0].E(new double[]{xhi - xlo, 0, 0});
                        edges[1].E(new double[]{xy, yhi - ylo, 0});
                        edges[2].E(new double[]{xz, yz, zhi - zlo});
                        System.out.println("edges: " + edges[0] + " " + edges[1] + " " + edges[2]);
                    }
                }
                item++;
            }
            else if (read == "atoms") {
                item++;
                String[] bits = line.split("[ \t]+");
                Vector p = box.getLeafList().getAtom(item - 1).getPosition();
                p.E(0);
                for (int i = 0; i < 3; i++) {
                    p.PEa1Tv1(Double.parseDouble(bits[i + 2]) - 0.5, edges[i]);
                }
                if (item == numAtoms) {
                    if (which == config) {
                        SimulationGraphic graphic = new SimulationGraphic(sim, space, sim.getController());
                        final DisplayBox display = new DisplayBox(sim, box, space, sim.getController());
                        graphic.add(display);

                        graphic.makeAndDisplayFrame();
                        return;
                    }
                    config++;
                }
            }
        }
    }

    public static class ConfigFromFileLAMMPSParam extends ParameterBase {
        public String filename;
        public int configNum = 0;
    }
}
