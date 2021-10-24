/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;


import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.graphics.ColorScheme;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceCheckBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorBoxFasterer;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.modifier.ModifierBoolean;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.simulation.Simulation;
import etomica.space.*;
import etomica.space3d.Tensor3D;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import javax.swing.*;
import java.awt.*;


public class StrainGraphic extends Simulation {

    public final CoordinateDefinitionLeaf coordinateDefinition;

    public Box box;
    public Boundary boundary;
    public Basis basis;
    public Primitive primitive;
    public SpeciesGeneral species;
    public IntegratorStrain integrator;

    public StrainGraphic(Space _space, int numAtoms, double density) {
        super(_space);
        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species);

        // TARGET
        double L = Math.pow(4.0 / density, 1.0 / 3.0);
        int n = (int) Math.round(Math.pow(numAtoms / 4, 1.0 / 3.0));
        boundary = new BoundaryDeformablePeriodic(space, n * L);
        box = this.makeBox(boundary);
        box.setNMolecules(species, numAtoms);

        primitive = new PrimitiveCubic(space, n * L);

        int[] nCells = new int[]{n, n, n};
        Basis basisFCC = new BasisCubicFcc();
        basis = new BasisBigCell(space, basisFCC, nCells);

        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(new int[]{1, 1, 1});

        integrator = new IntegratorStrain(box, n*L, 0.005);
        getController().addActivity(new ActivityIntegrate(integrator), Integer.MAX_VALUE, 1);
    }

    /**
     * @param args filename containing simulation parameters
     */
    public static void main(String[] args) {
        //set up simulation parameters
        SimOverlapParam params = new SimOverlapParam();
        if (args.length == 0) {
            params.numAtoms = 500;
            params.density = 1;
        }
        else {
            ParseArgs.doParseArgs(params, args);
        }
        double density = params.density;
        final int numAtoms = params.numAtoms;

        //instantiate simulation
        final StrainGraphic sim = new StrainGraphic(Space.getInstance(3), numAtoms, density);
        SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);

        JPanel cbPanel = new JPanel(new GridLayout(3,3));

        DeviceCheckBox cb11 = new DeviceCheckBox(sim.getController(), "xx", new ModifierStrain(1,1, sim.integrator.strain));
        cbPanel.add(cb11.graphic());
        DeviceCheckBox cb12 = new DeviceCheckBox(sim.getController(), "xy", new ModifierStrain(1,2, sim.integrator.strain));
        cbPanel.add(cb12.graphic());
        DeviceCheckBox cb13 = new DeviceCheckBox(sim.getController(), "xz", new ModifierStrain(1,3, sim.integrator.strain));
        cbPanel.add(cb13.graphic());

        DeviceCheckBox cb21 = new DeviceCheckBox(sim.getController(), "yx", new ModifierStrain(2,1, sim.integrator.strain));
        cbPanel.add(cb21.graphic());
        DeviceCheckBox cb22 = new DeviceCheckBox(sim.getController(), "yy", new ModifierStrain(2,2, sim.integrator.strain));
        cbPanel.add(cb22.graphic());
        DeviceCheckBox cb23 = new DeviceCheckBox(sim.getController(), "yz", new ModifierStrain(2,3, sim.integrator.strain));
        cbPanel.add(cb23.graphic());

        DeviceCheckBox cb31 = new DeviceCheckBox(sim.getController(), "zx", new ModifierStrain(3,1, sim.integrator.strain));
        cbPanel.add(cb31.graphic());
        DeviceCheckBox cb32 = new DeviceCheckBox(sim.getController(), "zy", new ModifierStrain(3,2, sim.integrator.strain));
        cbPanel.add(cb32.graphic());
        DeviceCheckBox cb33 = new DeviceCheckBox(sim.getController(), "zz", new ModifierStrain(3,3, sim.integrator.strain));
        cbPanel.add(cb33.graphic());

        simGraphic.getPanel().controlPanel.add(cbPanel);

        ColorSchemeByType csWhite = (ColorSchemeByType) simGraphic.getDisplayBox(sim.box).getColorScheme();
        csWhite.setColor(sim.species.getLeafType(), Color.WHITE);
        ColorSchemeStrain csStrain = new ColorSchemeStrain(sim.coordinateDefinition, 5);
        simGraphic.getDisplayBox(sim.box).setColorScheme(csStrain);

        DeviceCheckBox cbColor = new DeviceCheckBox(sim.getController(), "color", new ModifierBoolean() {
            boolean colorOn = true;

            @Override
            public void setBoolean(boolean b) {
                colorOn = b;
                simGraphic.getDisplayBox(sim.box).setColorScheme(colorOn ? csStrain : csWhite);
            }

            @Override
            public boolean getBoolean() {
                return colorOn;
            }
        });
        simGraphic.add(cbColor);

        simGraphic.makeAndDisplayFrame("Strain");
    }

    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int numAtoms = 256;
        public double density = 1;
    }

    public static class ColorSchemeStrain extends ColorScheme {

        public final CoordinateDefinition coordinateDefinition;
        public final double fac;

        public ColorSchemeStrain(CoordinateDefinition coordinateDefinition, double f) {
            this.coordinateDefinition = coordinateDefinition;
            fac = f;
        }

        @Override
        public Color getAtomColor(IAtom a) {
            Vector dr = coordinateDefinition.space.makeVector();
            dr.Ev1Mv2(a.getPosition(), coordinateDefinition.getLatticePosition(a));
            float r = Math.round(10*Math.max(0.0f, 1.0 - Math.abs(dr.getX(0)*fac)))*0.1f;
            float g = Math.round(10*Math.max(0.0f, 1.0 - Math.abs(dr.getX(1)*fac)))*0.1f;
            float b = Math.round(10*Math.max(0.0f, 1.0 - Math.abs(dr.getX(2)*fac)))*0.1f;
            return new Color(r,g,b);
        }
    }

    public static class ModifierStrain implements ModifierBoolean {

        public final int i, j;
        public final boolean[][] strain;

        public ModifierStrain(int i, int j, boolean[][] strain) {
            this.i = i - 1;
            this.j = j - 1;
            this.strain = strain;
        }

        @Override
        public void setBoolean(boolean b) {
            strain[i][j] = b;
        }

        @Override
        public boolean getBoolean() {
            return strain[i][j];
        }
    }

    public class IntegratorStrain extends IntegratorBoxFasterer {

        public final double L;
        public final boolean[][] strain;
        public double t;
        public double dt;
        public double amplitude;

        public IntegratorStrain(Box box, double L, double dt) {
            super(new PotentialComputeAggregate(), 1.0, box);
            this.L = L;
            this.strain = new boolean[3][3];
            this.dt = dt;
            this.amplitude = 0.05;
        }

        @Override
        protected void doStepInternal() {
            Tensor h = new Tensor3D(new double[][]{{L,0,0},{0,L,0},{0,0,L}});
            Tensor A = new Tensor3D();
            double y = Math.sin(t) * amplitude;
            for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
                    double x = i==j ? 1 : 0;
                    if (strain[i][j]) x += y;
                    A.setComponent(i, j, x);
                }
            }
            h.TE(A);
            ((BoundaryDeformablePeriodic)box.getBoundary()).setDimensions(h);
            for (IAtom a : box.getLeafList()) {
                Vector r = a.getPosition();
                r.E(coordinateDefinition.getLatticePosition(a));
                A.transform(r);
            }
            t += dt;
        }
    }
}
