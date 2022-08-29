/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.normalmode;

import etomica.action.IAction;
import etomica.atom.DiameterHash;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.graphics.ColorScheme;
import etomica.graphics.DeviceButton;
import etomica.graphics.DeviceSlider;
import etomica.graphics.SimulationGraphic;
import etomica.modifier.Modifier;
import etomica.space.Vector;
import etomica.space2d.Space2D;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Null;
import etomica.units.dimensions.Temperature;
import etomica.util.ParseArgs;

import java.awt.*;
import java.util.Arrays;

public class PathIntegralHOGraphic {

    public static void main(String[] args) {
        SimQuantumAO.OctaneParams params = new SimQuantumAO.OctaneParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            params.nBeads = 128;
        }

        double temperature = params.temperature;
        double k2 = params.k2;
        double k4 = params.k4;
        int nBeads = params.nBeads;
        boolean isTIA = params.isTIA;

        double omegaN = nBeads*temperature/SimQuantumAO.hbar; // 1/hbar*betan

        double massssss = 1.0;
        double omega2 = k2/massssss;
        if (isTIA){
            omega2 = omega2*(1.0 + omega2/12.0/omegaN/omegaN);
        }
        double beta = 1 / temperature;
        double betaN = beta/nBeads;

        final SimQuantumAO sim = new SimQuantumAO(Space2D.getInstance(), nBeads, temperature, k2, k4, omega2, isTIA, true);

        SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);

        double[] lambdaN = new double[nBeads];
        double[][] eigenvectors = new double[nBeads][nBeads];
        double[] q = new double[2*nBeads];  // initial q should all be 0
        for (int k=0; k<nBeads; k++) {
            double s = Math.sin(Math.PI*k/nBeads);
            lambdaN[k] = betaN * massssss*(4.0*omegaN*omegaN*s*s + omega2);
            for (int i = 0; i < nBeads; i++) {
                double arg = 2.0*Math.PI/nBeads*i*k;
                double cs = k > nBeads / 2 ? Math.sin(arg) : Math.cos(arg);
                double onetwo = k==0 || 2*k == nBeads ? 1 : 2;
                eigenvectors[i][k] = cs * Math.sqrt(onetwo / nBeads);
            }
        }

        ColorScheme colorScheme = new ColorScheme() {
            protected Color[] allColors;

            public Color getAtomColor(IAtom a) {
                int n = sim.box.getLeafList().size();
                if (allColors==null) {
                    allColors = new Color[768];
                    for (int i=0; i<256; i++) {
                        allColors[i] = new Color(255-i,i,0);
                    }
                    for (int i=0; i<256; i++) {
                        allColors[i+256] = new Color(0,255-i,i);
                    }
                    for (int i=0; i<256; i++) {
                        allColors[i+512] = new Color(i,0,255-i);
                    }
                }
                return allColors[(768*a.getLeafIndex()/n)];
            }
        };
        simGraphic.getDisplayBox(sim.box).setColorScheme(colorScheme);
        simGraphic.remove(simGraphic.getTrioControllerButton());
        DeviceButton newConfigButton = new DeviceButton(sim.getController(), new IAction() {
            @Override
            public void actionPerformed() {
                sim.atomMove.doTrial();
                sim.atomMove.acceptNotify();
                updateQ(sim.box.getLeafList(), q, eigenvectors);
                simGraphic.getDisplayBox(sim.box).repaint();
            }
        });
        newConfigButton.setLabel("New Config");
        simGraphic.getPanel().graphicsPanel.add(newConfigButton.graphic(), BorderLayout.SOUTH);
        DiameterHash diameterHash = simGraphic.getDisplayBox(sim.box).getDiameterHash();
        ((DiameterHashByType)diameterHash).setDiameter(sim.species().getLeafType(), 0.1);

        IAction updateLambda = new IAction() {
            @Override
            public void actionPerformed() {
                updateLambda(sim.box(), q, lambdaN, eigenvectors, sim.p1ah.getSpringConstants()[0], sim.integrator.getTemperature());
                simGraphic.getDisplayBox(sim.box).repaint();
            }
        };

        DeviceSlider k2Slider = new DeviceSlider(sim.getController(), new ModifierK2(massssss, isTIA, omegaN, sim));
        k2Slider.setMinimum(-2);
        k2Slider.setMaximum(4);
        k2Slider.setPrecision(1);
        k2Slider.setNMajor(6);
        k2Slider.setLabel("log10(k2)");
        k2Slider.setShowBorder(true);
        k2Slider.setPostAction(updateLambda);
        simGraphic.add(k2Slider);

        DeviceSlider temperatureSlider = new DeviceSlider(sim.getController(), new ModifierTemperature(nBeads, sim));
        temperatureSlider.setPrecision(1);
        temperatureSlider.setNMajor(5);
        temperatureSlider.setMaximum(2);
        temperatureSlider.setLabel("Temperature");
        temperatureSlider.setShowBorder(true);
        temperatureSlider.setPostAction(updateLambda);
        simGraphic.add(temperatureSlider);

        simGraphic.makeAndDisplayFrame("Path Integral Demo");

    }

    public static void updateLambda(Box box, double[] q, double[] lambda, double[][] eigenvectors, double k2, double temperature) {
        IAtomList atoms = box.getLeafList();
        int n = lambda.length;
        double beta = 1.0/(temperature);
        double betaN = beta/n;
        double omegaN = 1.0/(SimQuantumAO.hbar*betaN);
        double omega2 = k2/1;
        for (int k=0; k<n; k++) {
            double s = Math.sin(Math.PI * k / n);
            double newLambda = betaN * 1 * (4.0 * omegaN * omegaN * s * s + omega2);
            for (int d=0; d<q.length/n; d++) {
                q[d*n+k] *= Math.sqrt(lambda[k]/newLambda);
            }
            lambda[k] = newLambda;
        }

        updateR(atoms, q, eigenvectors);
    }

    public static void updateQ(IAtomList atoms, double[] q, double[][] eigenvectors) {
        Arrays.fill(q, 0);
        int n = eigenvectors.length;
        for (IAtom a : atoms) {
            int i = a.getLeafIndex();
            Vector r = a.getPosition();
            for (int d=0; d<r.getD(); d++) {
                for (int k=d*n; k<(d+1)*n; k++) {
                    q[k] += r.getX(d)*eigenvectors[i][k%n];
                }
            }
        }
    }

    public static void updateR(IAtomList atoms, double[] q, double[][] eigenvectors) {
        int n = eigenvectors.length;
        for (IAtom a : atoms) {
            int i = a.getLeafIndex();
            Vector r = a.getPosition();
            for (int d=0; d<r.getD(); d++) {
                double x = 0;
                for (int k=d*n; k<(d+1)*n; k++) {
                    x += eigenvectors[i][k%n]*q[k];
                }
                r.setX(d, x);
            }
        }
    }


    private static class ModifierTemperature implements Modifier {
        private final int nBeads;
        private final SimQuantumAO sim;

        public ModifierTemperature(int nBeads, SimQuantumAO sim) {
            this.nBeads = nBeads;
            this.sim = sim;
        }

        @Override
        public void setValue(double newValue) {
            if (newValue == 0) throw new IllegalArgumentException("T can't be 0");
            double beta = 1.0/(newValue);
            double betaN = beta/ nBeads;
            double omegaN = 1.0/(SimQuantumAO.hbar*betaN);

            sim.k2_kin = nBeads == 1 ? 0 : (1*omegaN*omegaN/ nBeads);
            // p2Bond.setk2
            ((MCMoveHOReal) sim.atomMove).setTemperature(newValue);
            // doesn't really matter, but we just use this to hold T
            sim.integrator.setTemperature(newValue);
        }

        @Override
        public double getValue() {
            return sim.integrator.getTemperature();
        }

        @Override
        public Dimension getDimension() {
            return Temperature.DIMENSION;
        }

        @Override
        public String getLabel() {
            return "Temperature";
        }
    }

    private static class ModifierK2 implements Modifier {
        private final double massssss;
        private final boolean isTIA;
        private final double omegaN;
        private final SimQuantumAO sim;

        public ModifierK2(double massssss, boolean isTIA, double omegaN, SimQuantumAO sim) {
            this.massssss = massssss;
            this.isTIA = isTIA;
            this.omegaN = omegaN;
            this.sim = sim;
        }

        @Override
        public void setValue(double newValue) {
            double k2new = Math.exp(newValue);
            double omega2 = k2new/ massssss;
            if (isTIA){
                omega2 = omega2*(1.0 + omega2/12.0/ omegaN / omegaN);
            }
            ((MCMoveHOReal) sim.atomMove).setOmega2(omega2);

            double[] k24 = sim.p1ah.getSpringConstants();
            sim.p1ah.setSpringConstants(k2new, k24[1]);
        }

        @Override
        public double getValue() {
            return Math.log(sim.p1ah.getSpringConstants()[0]);
        }

        @Override
        public Dimension getDimension() {
            return Null.DIMENSION;
        }

        @Override
        public String getLabel() {
            return "k2";
        }
    }
}
