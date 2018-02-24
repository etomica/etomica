/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.freeenergy.npath;

import etomica.action.IAction;
import etomica.atom.*;
import etomica.box.Box;
import etomica.data.DataSourceIndependentSimple;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataSink;
import etomica.data.meter.MeterStructureFactor;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.graphics.*;
import etomica.math.numerical.ArrayReader1D;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierBoolean;
import etomica.simulation.Simulation;
import etomica.space.*;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresRotating;
import etomica.units.dimensions.*;
import etomica.units.dimensions.Dimension;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;
import java.io.*;
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
        boolean doSfac = params.doSfac;
        double cutoffS = params.cutS;
        double thresholdS = params.thresholdS;
        boolean GUI = params.GUI;
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
        Species species = new SpeciesSpheresRotating(sim, space);
        sim.addSpecies(species);
        Box box = null;
        int config = 0;
        List<Vector[]> allCoords = new ArrayList<>();
        Vector[] iCoords = null;
        List<Vector[]> allEdges = new ArrayList<>();
        List<DataFunction> sfacData = new ArrayList<>();
        List<DataFunction.DataInfoFunction> sfacDataInfo = new ArrayList<>();
        DataTag sfacTag = new DataTag();
        ModifierConfiguration modifierConfig =  null;
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
                    if (config == 0) {
                        Boundary boundary = new BoundaryDeformablePeriodic(space, edges);
                        box = sim.makeBox(boundary);
                        box.setNMolecules(species, numAtoms);
                        modifierConfig = new ModifierConfiguration(box, allCoords, allEdges, sfacData, sfacDataInfo, sfacTag);
                    }
                    iCoords = new Vector[numAtoms];
                }
                item = 0;
                continue;
            }
            if (read.equals("numAtoms")) {
                numAtoms = Integer.parseInt(line);
                System.out.println("numAtoms: " + numAtoms);
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
                    if (doSfac) {
                        File sfile = new File(config + ".sfac");
                        if (!sfile.exists()) {
                            modifierConfig.setValue(config);
                            MeterStructureFactor meter = new MeterStructureFactor(space, box, cutoffS);
                            writeFile(meter, thresholdS, config);
                        }
                    }
                    config++;
                    edges = space.makeVectorArray(3);

                }
            }
        }
        if (!GUI) return;
        SimulationGraphic graphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
        final DisplayBox display = new DisplayBox(sim, box);
        final DisplayBoxCanvasG3DSys.OrientedSite site = new DisplayBoxCanvasG3DSys.OrientedSite(0.5, Color.WHITE, 0.5);
        ((DisplayBoxCanvasG3DSys) display.canvas).setOrientationSites((AtomTypeOriented) species.getAtomType(0), new DisplayBoxCanvasG3DSys.OrientedSite[]{site});
        graphic.add(display);

        ColorSchemeDeviation colorScheme = new ColorSchemeDeviation(box, space);
        graphic.getDisplayBox(box).setColorScheme(colorScheme);
        final AtomFilterDeviation atomFilter = new AtomFilterDeviation(colorScheme);
        display.setAtomFilter(atomFilter);
        modifierConfig.setColorScheme(colorScheme);

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
        filterSlider.setMaximum(2);
        filterSlider.setPrecision(2);
        filterSlider.setNMajor(4);
        Box finalBox = box;
        filterSlider.setPostAction(() -> graphic.getDisplayBox(finalBox).repaint());
        graphic.add(filterSlider);
        DeviceSlider configSlider = new DeviceSlider(sim.getController(), modifierConfig);
        configSlider.setMaximum(config - 1);
        configSlider.setNMajor(config);
        configSlider.setMinimum(0);
        configSlider.setPrecision(0);
        Box finalBox1 = box;
        modifierConfig.setPostAction(new IAction() {

            @Override
            public void actionPerformed() {
                graphic.getDisplayBox(finalBox1).repaint();
            }
        });
        graphic.add(configSlider);
        configSlider.setValue(0);

        DisplayPlot sfacPlot = null;
        IDataSink sfacPlotSink = null;
        for (int i = 0; i <= config; i++) {
            int xMin = (int) cutoffS + 1;
            File sfile = new File(i + ".sfac");
            if (sfile.exists()) {
                double[][] xy = ArrayReader1D.getFromFile(i + ".sfac");
                double[] x = new double[xy.length];
                double[] y = new double[xy.length];
                for (int j = 0; j < x.length; j++) {
                    x[j] = xy[j][0];
                    if (xMin > x[j]) xMin = (int) x[j];
                    y[j] = xy[j][1];
                }
                DataFunction data = new DataFunction(new int[]{y.length}, y);
                DataDoubleArray.DataInfoDoubleArray xDataInfo = new DataDoubleArray.DataInfoDoubleArray("wave vector", new CompoundDimension(new Dimension[]{Length.DIMENSION}, new double[]{-1}), new int[]{x.length});
                DataSourceIndependentSimple xDataSource = new DataSourceIndependentSimple(x, xDataInfo);
                DataFunction.DataInfoFunction dataInfo = new DataFunction.DataInfoFunction("structure factor", Null.DIMENSION, xDataSource);
                dataInfo.addTag(sfacTag);

                if (sfacPlot == null) {
                    sfacPlot = new DisplayPlot();
                    sfacPlot.getPlot().setYLog(true);
                    sfacPlot.getPlot().setXRange(xMin, cutoffS);
                    sfacPlot.getPlot().setYRange(Math.log10(thresholdS), Math.log10(2));
                    sfacPlot.setLabel("structure factor");
                    sfacPlot.setDoLegend(false);
                    sfacPlotSink = sfacPlot.getDataSet().makeDataSink();
                    sfacPlotSink.putDataInfo(dataInfo);
                    sfacPlot.setDoDrawLines(new DataTag[]{sfacTag}, false);
                    graphic.add(sfacPlot);
                    modifierConfig.setSfacPlotSink(sfacPlotSink);
                }
                sfacData.add(data);
                sfacDataInfo.add(dataInfo);
            }
        }

        modifierConfig.setValue(0);
        ModifierConfiguration finalModifierConfig = modifierConfig;
        DeviceToggleRadioButtons deviceDoPrevious = new DeviceToggleRadioButtons(new ModifierBoolean() {
            @Override
            public void setBoolean(boolean b) {
                finalModifierConfig.setDoPrevious(b);
            }
    
            @Override
            public boolean getBoolean() {
                return finalModifierConfig.getDoPrevious();
            }
        },"Deviation from","previous","original");
        graphic.add(deviceDoPrevious);

        DeviceCheckBox deviceShowDirection = new DeviceCheckBox("Show direction", new ModifierBoolean() {
            boolean sitesShown = true;

            @Override
            public boolean getBoolean() {
                return sitesShown;
            }

            @Override
            public void setBoolean(boolean b) {
                if (b == sitesShown) return;
                if (b) {
                    ((DisplayBoxCanvasG3DSys) display.canvas).setOrientationSites((AtomTypeOriented) species.getAtomType(0), new DisplayBoxCanvasG3DSys.OrientedSite[]{site});
                } else {
                    ((DisplayBoxCanvasG3DSys) display.canvas).setOrientationSites((AtomTypeOriented) species.getAtomType(0), new DisplayBoxCanvasG3DSys.OrientedSite[0]);
                }
                sitesShown = b;
                display.repaint();
            }
        });
        graphic.add(deviceShowDirection);

        graphic.makeAndDisplayFrame();
    }

    public static void writeFile(MeterStructureFactor meter, double threshold, int config) throws IOException {
        IData data = meter.getData();
        IData xData = ((DataFunction.DataInfoFunction) meter.getDataInfo()).getXDataSource().getIndependentData(0);
        FileWriter fw = new FileWriter(config + ".sfac");
        for (int i = 0; i < data.getLength(); i++) {
            double sfac = data.getValue(i);
            if (sfac > threshold) fw.write(xData.getValue(i) + " " + sfac + "\n");
        }
        fw.close();
    }

    public static class ConfigFromFileLAMMPSParam extends ParameterBase {
        public String filename;
        public boolean doSfac = false;
        public double cutS = 8;
        public double thresholdS = 0.001;
        public boolean GUI = true;
    }

    public static class ColorSchemeDeviation extends ColorScheme implements ColorSchemeCollective {
        protected final Vector[] scaledCoords0, rescaledCoords0;
        protected final Box box;
        protected final Vector dr;
        protected final Color[] colors;
        protected final Vector[] edges;
        protected final Tensor t;
        protected final Vector r0;
        protected double rNbr;

        public ColorSchemeDeviation(Box box, Space space) {
            this.box = box;
            dr = space.makeVector();
            colors = new Color[511];
            for (int i = 0; i < 256; i++) {
                colors[i] = new Color(0, i, 255 - i);
            }
            for (int i = 1; i < 256; i++) {
                colors[255 + i] = new Color(i, 255 - i, 0);
            }
            edges = space.makeVectorArray(3);
            for (int i = 0; i < 3; i++) {
                edges[i] = box.getBoundary().getEdgeVector(i);
            }
            int numAtoms = box.getLeafList().size();
            t = space.makeTensor();
            scaledCoords0 = space.makeVectorArray(numAtoms);
            rescaledCoords0 = space.makeVectorArray(numAtoms);
            r0 = space.makeVector();
        }
        
        public void setLattice(Vector[] latticeEdges, Vector[] latticeCoords) {
            t.E(latticeEdges);
            t.invert();
            for (int i = 0; i < scaledCoords0.length; i++) {
                scaledCoords0[i].E(latticeCoords[i]);
                t.transform(scaledCoords0[i]);
            }
        }

        @Override
        public void colorAllAtoms() {
            IAtomList atoms = box.getLeafList();

            double vol = box.getBoundary().volume() / atoms.size();
            double a = Math.cbrt(vol);
            rNbr = a * Math.sqrt(3) / 2;

            BoundaryDeformablePeriodic boundary = (BoundaryDeformablePeriodic) box.getBoundary();
            for (int i = 0; i < 3; i++) {
                edges[i] = boundary.getEdgeVector(i);
            }
            t.E(edges);

            for (int i = 0; i < atoms.size(); i++) {
                r0.E(scaledCoords0[i]);
                t.transform(r0);
                rescaledCoords0[i].E(r0);
                dr.Ev1Mv2(atoms.get(i).getPosition(), r0);
                boundary.nearestImage(dr);
                double r2 = dr.squared();
            }

        }

        public Vector getDisplacement(IAtom a) {
            dr.Ev1Mv2(a.getPosition(), rescaledCoords0[a.getLeafIndex()]);
            box.getBoundary().nearestImage(dr);
            return dr;
        }

        public double getRelativeDisplacement(IAtom a) {
            return Math.sqrt(getDisplacement(a).squared()) / rNbr;
        }

        @Override
        public Color getAtomColor(IAtom a) {
            double s = getRelativeDisplacement(a);
            if (s >= 2) return colors[510];
            int i = (int) (510.9999 * s / 2);
            return colors[i];
        }
    }

    public static class ModifierConfiguration implements Modifier {

        protected final List<Vector[]> allCoords;
        protected final List<Vector[]> allEdges;
        protected final List<DataFunction> sfacData;
        protected final List<DataFunction.DataInfoFunction> sfacDataInfo;
        protected final DataFunction.DataInfoFunction emptySfacDataInfo;
        protected final DataFunction emptySfacData;
        protected int configIndex = -1;
        protected Box box;
        protected IDataSink sfacPlotSink;
        protected ColorSchemeDeviation colorScheme;
        protected boolean doPrevious;
        protected IAction postAction;

        public ModifierConfiguration(Box box, List<Vector[]> allCoords, List<Vector[]> allEdges, List<DataFunction> sfacData, List<DataFunction.DataInfoFunction> sfacDataInfo, DataTag sfacTag) {
            this.box = box;
            this.allCoords = allCoords;
            this.allEdges = allEdges;
            this.sfacData = sfacData;
            this.sfacDataInfo = sfacDataInfo;
            this.emptySfacData = new DataFunction(new int[]{0}, new double[0]);
            DataDoubleArray.DataInfoDoubleArray xDataInfo = new DataDoubleArray.DataInfoDoubleArray("wave vector", new CompoundDimension(new Dimension[]{Length.DIMENSION}, new double[]{-1}), new int[]{0});
            double[] x = new double[0];
            DataSourceIndependentSimple xDataSource = new DataSourceIndependentSimple(x, xDataInfo);
            emptySfacDataInfo = new DataFunction.DataInfoFunction("structure factor", Null.DIMENSION, xDataSource);
            emptySfacDataInfo.addTag(sfacTag);
        }
        
        public void setPostAction(IAction postAction) {
            this.postAction = postAction;
        }

        public boolean getDoPrevious() {
            return doPrevious;
        }
        
        public void setDoPrevious(boolean doPrevious) {
            if (this.doPrevious == doPrevious) return;
            this.doPrevious = doPrevious;
            int foo = configIndex;
            configIndex = -1;
            setValue(foo);
        }

        public void setColorScheme(ColorSchemeDeviation colorScheme) {
            this.colorScheme = colorScheme;
        }

        public void setSfacPlotSink(IDataSink sfacPlotSink) {
            this.sfacPlotSink = sfacPlotSink;
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
            for (int i = 0; i < atoms.size(); i++) {
                IAtom a = atoms.get(i);
                a.getPosition().E(myCoords[i]);
                if (colorScheme != null) {
                    Vector orientation = colorScheme.getDisplacement(a);
                    orientation.normalize();
                    ((IAtomOriented) a).getOrientation().setDirection(orientation);
                }
            }
            if (colorScheme != null) {
                if (doPrevious && configIndex > 0) {
                    colorScheme.setLattice(allEdges.get(configIndex - 1), allCoords.get(configIndex - 1));
                } else {
                    colorScheme.setLattice(allEdges.get(0), allCoords.get(0));
                }
            }

            if (sfacPlotSink != null) {
                if (configIndex < sfacData.size()) {
                    sfacPlotSink.putDataInfo(sfacDataInfo.get(configIndex));
                    sfacPlotSink.putData(sfacData.get(configIndex));
                } else {
                    sfacPlotSink.putDataInfo(emptySfacDataInfo);
                    sfacPlotSink.putData(emptySfacData);
                }
            }
            if (postAction!= null) {
                postAction.actionPerformed();
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
        public boolean test(IAtom a) {
            double x = colorScheme.getRelativeDisplacement(a);
            return x >= threshold;
        }
    }

}
