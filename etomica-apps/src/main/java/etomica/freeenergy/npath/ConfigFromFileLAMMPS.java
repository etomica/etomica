package etomica.freeenergy.npath;

import etomica.action.IAction;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.graphics.*;
import etomica.modifier.Modifier;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Dimension;
import etomica.units.Quantity;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * This reads ATOM files from LAMMPS and creates a display of one of the
 * configurations.
 * <p>
 * Created by andrew on 4/25/17.
 */
public class ConfigFromFileLAMMPS {

    public static void main(String[] args) throws IOException {
        ConfigFromFileLAMMPSParam params = new ConfigFromFileLAMMPSParam();
        if (args.length == 0) {
            throw new RuntimeException("Usage: ConfigFromFileLAMMPS -filename sim.atom -configNum n -crystal CRYSTAL");
        }
        ParseArgs.doParseArgs(params, args);
        String filename = params.filename;
        boolean colorDeviation = params.colorDeviation;
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
        List<Vector[]> allCoords = new ArrayList<>();
        Vector[] iCoords = null;
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
                    iCoords = new Vector[numAtoms];
                }
                item = 0;
                continue;
            }
            if (read.equals("numAtoms")) {
                numAtoms = Integer.parseInt(line);
                System.out.println("numAtoms: " + numAtoms);
                box.setNMolecules(species, numAtoms);
            } else if (read.equals("boundary")) {
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
            } else if (read.equals("atoms")) {
                item++;
                String[] bits = line.split("[ \t]+");
                Vector p = space.makeVector();
                p.E(0);
                for (int i = 0; i < 3; i++) {
                    p.PEa1Tv1(Double.parseDouble(bits[i + 2]) - 0.5, edges[i]);
                }
                iCoords[item - 1] = p;
                if (item == numAtoms) {
                    allCoords.add(iCoords);
                    config++;
                }
            }
        }
        SimulationGraphic graphic = new SimulationGraphic(sim, space, sim.getController());
        final DisplayBox display = new DisplayBox(sim, box, space, sim.getController());
        graphic.add(display);
        if (colorDeviation) {
            graphic.getDisplayBox(box).setColorScheme(new ColorSchemeDeviation(allCoords.get(0), box, space));
        }
        ModifierConfiguration modifierConfig = new ModifierConfiguration(box, allCoords);
        DeviceSlider configSlider = new DeviceSlider(sim.getController(), modifierConfig);
        configSlider.setMaximum(config - 2);
        configSlider.setNMajor(config - 1);
        configSlider.setMinimum(0);
        configSlider.setPrecision(0);
        configSlider.setPostAction(new IAction() {

            @Override
            public void actionPerformed() {
                graphic.getDisplayBox(box).repaint();
            }
        });
        graphic.add(configSlider);

        graphic.makeAndDisplayFrame();
        modifierConfig.setValue(0);
    }

    public static class ConfigFromFileLAMMPSParam extends ParameterBase {
        public String filename;
        public boolean colorDeviation = true;
    }

    public static class ColorSchemeDeviation extends ColorScheme implements ColorSchemeCollective {
        protected final Vector[] latticeCoords;
        protected final Box box;
        protected final Vector dr;
        protected final Color[] colors;
        protected double rMax;

        public ColorSchemeDeviation(Vector[] latticeCoords, Box box, Space space) {
            this.latticeCoords = latticeCoords;
            this.box = box;
            dr = space.makeVector();
            colors = new Color[256];
            for (int i = 0; i < 256; i++) {
                colors[i] = new Color(i, 0, 255 - i);
            }
        }

        @Override
        public void colorAllAtoms() {
            double maxR2 = 0;
            IAtomList atoms = box.getLeafList();
            Boundary boundary = box.getBoundary();
            for (int i = 0; i < atoms.getAtomCount(); i++) {
                dr.Ev1Mv2(atoms.getAtom(i).getPosition(), latticeCoords[i]);
                boundary.nearestImage(dr);
                double r2 = dr.squared();
                if (r2 > maxR2) maxR2 = r2;
            }
            rMax = maxR2 == 0 ? 1 : Math.sqrt(maxR2);
        }

        @Override
        public Color getAtomColor(IAtom a) {
            dr.Ev1Mv2(a.getPosition(), latticeCoords[a.getLeafIndex()]);
            box.getBoundary().nearestImage(dr);
            int i = (int) (255.9999 * (Math.sqrt(dr.squared()) / rMax));
            return colors[i];
        }
    }

    public static class ModifierConfiguration implements Modifier {

        protected int configIndex = -1;
        protected Box box;
        protected List<Vector[]> allCoords;

        public ModifierConfiguration(Box box, List<Vector[]> allCoords) {
            this.box = box;
            this.allCoords = allCoords;
        }

        @Override
        public double getValue() {
            return configIndex;
        }

        @Override
        public void setValue(double newValue) {
            int nv = (int) Math.round(newValue);
            if (nv == configIndex) return;
            configIndex = nv;
            IAtomList atoms = box.getLeafList();
            Vector[] myCoords = allCoords.get(configIndex);
            for (int i = 0; i < atoms.getAtomCount(); i++) {
                atoms.getAtom(i).getPosition().E(myCoords[i]);
            }
        }

        @Override
        public Dimension getDimension() {
            return Quantity.DIMENSION;
        }

        @Override
        public String getLabel() {
            return "Configuration";
        }
    }
}
