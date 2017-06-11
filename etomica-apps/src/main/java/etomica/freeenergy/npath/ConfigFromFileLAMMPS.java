package etomica.freeenergy.npath;

import etomica.action.IAction;
import etomica.atom.AtomFilter;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.graphics.*;
import etomica.modifier.Modifier;
import etomica.molecule.IMolecule;
import etomica.simulation.Simulation;
import etomica.space.*;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Dimension;
import etomica.units.Null;
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
        List<Vector[]> allEdges = new ArrayList<>();
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
                    allEdges.add(edges);
                    if (config == 1) {
                        Boundary boundary = new BoundaryDeformablePeriodic(space, edges);
                        box.setBoundary(boundary);
                    }
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
                    edges = space.makeVectorArray(3);
                }
            }
        }
        SimulationGraphic graphic = new SimulationGraphic(sim, space, sim.getController());
        final DisplayBox display = new DisplayBox(sim, box, space, sim.getController());
        graphic.add(display);
        if (colorDeviation) {
            ColorSchemeDeviation colorScheme = new ColorSchemeDeviation(allCoords.get(0), box, space);
            graphic.getDisplayBox(box).setColorScheme(colorScheme);
            final AtomFilterDeviation atomFilter = new AtomFilterDeviation(colorScheme);
            display.setAtomFilter(atomFilter);

            DeviceSlider filterSlider = new DeviceSlider(sim.getController(), new Modifier() {
                @Override
                public double getValue() {
                    return atomFilter.getThreshold();
                }

                @Override
                public void setValue(double newValue) {
                    atomFilter.setThreshold(newValue);
                }

                @Override
                public Dimension getDimension() {
                    return Null.DIMENSION;
                }

                @Override
                public String getLabel() {
                    return "Displacement Threshold";
                }
            });
            filterSlider.setMinimum(0);
            filterSlider.setMaximum(1);
            filterSlider.setPrecision(2);
            filterSlider.setNMajor(5);
            filterSlider.setPostAction(new IAction() {

                @Override
                public void actionPerformed() {
                    graphic.getDisplayBox(box).repaint();
                }
            });
            graphic.add(filterSlider);
        }
        ModifierConfiguration modifierConfig = new ModifierConfiguration(box, allCoords, allEdges);
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

        modifierConfig.setValue(0);

        graphic.makeAndDisplayFrame();
    }

    public static class ConfigFromFileLAMMPSParam extends ParameterBase {
        public String filename;
        public boolean colorDeviation = true;
    }

    public static class ColorSchemeDeviation extends ColorScheme implements ColorSchemeCollective {
        protected final Vector[] scaledCoords0, rescaledCoords0;
        protected final Box box;
        protected final Vector dr;
        protected final Color[] colors;
        protected final Vector[] edges0, edges;
        protected final Tensor t0, t;
        protected final Vector r0;
        protected double rMax;

        public ColorSchemeDeviation(Vector[] latticeCoords, Box box, Space space) {
            this.box = box;
            dr = space.makeVector();
            colors = new Color[256];
            for (int i = 0; i < 256; i++) {
                colors[i] = new Color(i, 0, 255 - i);
            }
            edges0 = space.makeVectorArray(3);
            for (int i = 0; i < 3; i++) {
                edges0[i] = box.getBoundary().getEdgeVector(i);
            }
            t0 = space.makeTensor();
            t0.E(edges0);
            t0.invert();
            edges = space.makeVectorArray(3);
            t = space.makeTensor();
            scaledCoords0 = space.makeVectorArray(latticeCoords.length);
            for (int i = 0; i < scaledCoords0.length; i++) {
                scaledCoords0[i].E(latticeCoords[i]);
                t0.transform(scaledCoords0[i]);
            }
            rescaledCoords0 = space.makeVectorArray(latticeCoords.length);
            r0 = space.makeVector();
        }

        @Override
        public void colorAllAtoms() {
            double maxR2 = 0;
            IAtomList atoms = box.getLeafList();
            BoundaryDeformablePeriodic boundary = (BoundaryDeformablePeriodic) box.getBoundary();
            for (int i = 0; i < 3; i++) {
                edges[i] = boundary.getEdgeVector(i);
            }
            t.E(edges);

            for (int i = 0; i < atoms.getAtomCount(); i++) {
                r0.E(scaledCoords0[i]);
                t.transform(r0);
                rescaledCoords0[i].E(r0);
                dr.Ev1Mv2(atoms.getAtom(i).getPosition(), r0);
                boundary.nearestImage(dr);
                double r2 = dr.squared();
                if (r2 > maxR2) maxR2 = r2;
            }
            rMax = maxR2 < 1e-16 ? 1 : Math.sqrt(maxR2);
        }

        public double getRelativeDisplacement(IAtom a) {
            dr.Ev1Mv2(a.getPosition(), rescaledCoords0[a.getLeafIndex()]);
            box.getBoundary().nearestImage(dr);
            return Math.sqrt(dr.squared()) / rMax;
        }

        @Override
        public Color getAtomColor(IAtom a) {
            int i = (int) (255.9999 * getRelativeDisplacement(a));
            return colors[i];
        }
    }

    public static class ModifierConfiguration implements Modifier {

        protected final List<Vector[]> allCoords;
        protected final List<Vector[]> allEdges;
        protected int configIndex = -1;
        protected Box box;

        public ModifierConfiguration(Box box, List<Vector[]> allCoords, List<Vector[]> allEdges) {
            this.box = box;
            this.allCoords = allCoords;
            this.allEdges = allEdges;
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
            for (int i = 0; i < 3; i++) {
                ((BoundaryDeformablePeriodic) box.getBoundary()).setEdgeVector(i, allEdges.get(configIndex)[i]);
            }
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

    public static class AtomFilterDeviation implements AtomFilter {

        protected final ColorSchemeDeviation colorScheme;
        protected double threshold = 0;

        public AtomFilterDeviation(ColorSchemeDeviation colorScheme) {
            this.colorScheme = colorScheme;
        }

        public double getThreshold() {
            return threshold;
        }

        public void setThreshold(double threshold) {
            this.threshold = threshold;
        }

        @Override
        public boolean accept(IAtom a) {
            double x = colorScheme.getRelativeDisplacement(a);
            return x >= threshold;
        }

        @Override
        public boolean accept(IMolecule mole) {
            return false;
        }
    }
}