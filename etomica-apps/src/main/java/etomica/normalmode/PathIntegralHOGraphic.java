/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.normalmode;

import etomica.action.IAction;
import etomica.atom.DiameterHash;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.graphics.*;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierBoolean;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space2d.Space2D;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Null;
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
        double hbar = params.hbar;
        double k2 = params.k2;
        double k4 = params.k4;
        int nBeads = params.nBeads;
        boolean isTIA = params.isTIA;

        double omegaN = Math.sqrt(nBeads)*temperature/hbar;

        double massssss = 1.0;
        double omega2 = k2/massssss;
        if (isTIA){
            omega2 = omega2*(1.0 + omega2/12.0/(nBeads*omegaN*omegaN));
        }

        final SimQuantumAO sim = new SimQuantumAO(Space2D.getInstance(), SimQuantumAO.MoveChoice.Real, 1.0, nBeads, temperature, k2, k4, omega2, isTIA, hbar);

        SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);

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
        TransformerNormalMode transformerNM = new TransformerNormalMode(sim.box, k2, temperature, hbar);
        TransformerReal transformerReal = new TransformerReal(sim.box, (MCMoveHOReal) sim.move);
        NewConfigAction newConfigAction = new NewConfigAction(sim, simGraphic, transformerNM);
        DeviceButton newConfigButton = new DeviceButton(sim.getController(), newConfigAction);
        newConfigButton.setLabel("New Config");
        simGraphic.getPanel().graphicsPanel.add(newConfigButton.graphic(), BorderLayout.SOUTH);
        DiameterHash diameterHash = simGraphic.getDisplayBox(sim.box).getDiameterHash();
        ((DiameterHashByType)diameterHash).setDiameter(sim.species().getLeafType(), 0.1);

        FullUpdateAction updateLambda = new FullUpdateAction(sim, simGraphic, transformerNM);

        DeviceSlider k2Slider = new DeviceSlider(sim.getController(), new ModifierK2(massssss, isTIA, omegaN, sim));
        k2Slider.setMinimum(-2);
        k2Slider.setMaximum(4);
        k2Slider.setPrecision(1);
        k2Slider.setNMajor(6);
        k2Slider.setLabel("log10(k2)");
        k2Slider.setShowBorder(true);
        k2Slider.doUpdate();
        k2Slider.setPostAction(updateLambda);
        simGraphic.add(k2Slider);

        DeviceSlider temperatureSlider = new DeviceSlider(sim.getController(), new ModifierTemperature(nBeads, sim, hbar));
        temperatureSlider.setMinimum(-2);
        temperatureSlider.setMaximum(2);
        temperatureSlider.setPrecision(1);
        temperatureSlider.setNMajor(4);
        temperatureSlider.setLabel("log10(T)");
        temperatureSlider.setShowBorder(true);
        temperatureSlider.doUpdate();
        temperatureSlider.setPostAction(updateLambda);
        simGraphic.add(temperatureSlider);

        DeviceCheckBox realCheckbox = new DeviceCheckBox(sim.getController(), "real space", new ModifierBoolean() {
            @Override
            public void setBoolean(boolean b) {
                Transformer t = b ? transformerReal : transformerNM;
                t.updateAndQ(sim.p1ah.getSpringConstants()[0], sim.integrator.getTemperature());
                newConfigAction.setTransformer(t);
                updateLambda.setTransformer(t);
            }

            @Override
            public boolean getBoolean() {
                return newConfigAction.transformer instanceof TransformerReal;
            }
        });
        simGraphic.add(realCheckbox);

        simGraphic.makeAndDisplayFrame("Path Integral Demo");

    }

    public static class TransformerReal implements Transformer {

        public final double[] q;
        public final MCMoveHOReal move;
        public final Box box;
        public final double[] chainSigmas;

        public TransformerReal(Box box, MCMoveHOReal move) {
            this.box = box;
            this.move = move;
            q = new double[box.getSpace().D()*box.getLeafList().size()];
            chainSigmas = move.getChainSigmas().clone();
            updateQ();
        }

        public void updateR() {
            int n = box.getLeafList().size();
            int D = box.getSpace().D();
            double[][] centerCoefficients = move.getCenterCoefficients();
            double[] R11 = centerCoefficients[0], R1N = centerCoefficients[1];
            IAtom atom0 = null, atomPrev = null;
            for (IAtom a : box.getLeafList()) {
                Vector ri = a.getPosition();
                int i = a.getLeafIndex();
                int N = n - i;
                if (i == 0) {
                    for (int d=0; d<D; d++) {
                        ri.setX(d, q[d*n]);
                    }
                    atom0 = atomPrev = a;
                    continue;
                }
                for (int d=0; d<D; d++) {
                    ri.setX(d, q[d*n+i]);
                }
                ri.PEa1Tv1(R11[N], atomPrev.getPosition());
                ri.PEa1Tv1(R1N[N], atom0.getPosition());
                atomPrev = a;
            }
        }

        @Override
        public void updateQ() {
            int n = box.getLeafList().size();
            int D = box.getSpace().D();
            double[][] centerCoefficients = move.getCenterCoefficients();
            double[] R11 = centerCoefficients[0], R1N = centerCoefficients[1];
            IAtom atom0 = null, atomPrev = null;
            for (IAtom a : box.getLeafList()) {
                int i = a.getLeafIndex();
                Vector r = a.getPosition();
                if (i == 0) {
                    for (int d=0; d<D; d++) {
                        q[d*n] = r.getX(d);
                    }
                    atom0 = atomPrev = a;
                    continue;
                }
                int N = n - i;
                Vector dr = box.getSpace().makeVector();;
                dr.E(a.getPosition());
                dr.PEa1Tv1(-R11[N], atomPrev.getPosition());
                dr.PEa1Tv1(-R1N[N], atom0.getPosition());
                for (int d=0; d<r.getD(); d++) {
                    q[d*n + i] = dr.getX(d);
                }
                atomPrev = a;
            }
        }

        @Override
        public void update(double k2, double temperature) {
            // do nothing.  we let MCMove get updated and then retrieve stuff from it.
            double[] newChainSigmas = move.getChainSigmas();
            int n = box.getLeafList().size();
            for (int d=0; d<box.getSpace().D(); d++) {
                for (int i=0; i<n; i++) {
                    if (i == 0) {
                        q[d*n] *= newChainSigmas[0] / chainSigmas[0];
                    }
                    else {
                        int N = n - i;
                        q[d*n + i] *= newChainSigmas[N] / chainSigmas[N];
                    }
                }
            }
            System.arraycopy(newChainSigmas, 0, chainSigmas, 0, chainSigmas.length);

            updateR();
        }

        public void updateAndQ(double k2, double temperature) {
            // do nothing.  we let MCMove get updated and then retrieve stuff from it.
            double[] newChainSigmas = move.getChainSigmas();
            System.arraycopy(newChainSigmas, 0, chainSigmas, 0, chainSigmas.length);

            updateQ();
        }

    }

    public static class TransformerNormalMode implements Transformer {

        public final double[] q;
        public final double[] lambda;
        public final double[][] eigenvectors;
        public final Box box;
        public final double hbar;

        public TransformerNormalMode(Box box, double k2, double temperature, double hbar) {
            this.box = box;
            this.hbar = hbar;
            int nBeads = box.getLeafList().size();
            double mass = 1;
            double omegaN = Math.sqrt(nBeads)*temperature/hbar;

            double omega2 = k2/mass;

            double beta = 1 / temperature;
            double betaN = beta/nBeads;
            lambda = new double[nBeads];
            eigenvectors = new double[nBeads][nBeads];
            q = new double[box.getSpace().D()*nBeads];  // initial q should all be 0
            for (int k=0; k<nBeads; k++) {
                double s = Math.sin(Math.PI*k/nBeads);
                lambda[k] = beta * mass * (4.0*omegaN*omegaN*s*s + omega2/nBeads);
                for (int i = 0; i < nBeads; i++) {
                    double arg = 2.0*Math.PI/nBeads*i*k;
                    double cs = k > nBeads / 2 ? Math.sin(arg) : Math.cos(arg);
                    double onetwo = k==0 || 2*k == nBeads ? 1 : 2;
                    eigenvectors[i][k] = cs * Math.sqrt(onetwo / nBeads);
                }
            }
        }

        public void updateR() {
            int n = eigenvectors.length;
            for (IAtom a : box.getLeafList()) {
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

        @Override
        public void updateQ() {
            Arrays.fill(q, 0);
            int n = eigenvectors.length;
            for (IAtom a : box.getLeafList()) {
                int i = a.getLeafIndex();
                Vector r = a.getPosition();
                for (int d=0; d<r.getD(); d++) {
                    for (int k=d*n; k<(d+1)*n; k++) {
                        q[k] += r.getX(d)*eigenvectors[i][k%n];
                    }
                }
            }
        }

        protected void updateInternal(double k2, double temperature) {
            int n = lambda.length;
            double beta = 1.0/(temperature);
            double betaN = beta/n;
            double omegaN = Math.sqrt(n)/(beta*hbar);
            double omega2 = k2;
            for (int k=0; k<n; k++) {
                double s = Math.sin(Math.PI * k / n);
                double newLambda = beta * (4.0 * omegaN * omegaN * s * s + omega2/n);
                for (int d=0; d<q.length/n; d++) {
                    q[d*n+k] *= Math.sqrt(lambda[k]/newLambda);
                }
                lambda[k] = newLambda;
            }
        }

        @Override
        public void update(double k2, double temperature) {
            int n = lambda.length;
            double beta = 1.0/(temperature);
            double betaN = beta/n;
            double omegaN = Math.sqrt(n)/(beta*hbar);
            double omega2 = k2;
            for (int k=0; k<n; k++) {
                double s = Math.sin(Math.PI * k / n);
                double newLambda = beta * 1 * (4.0 * omegaN * omegaN * s * s + omega2/n);
                for (int d=0; d<q.length/n; d++) {
                    q[d*n+k] *= Math.sqrt(lambda[k]/newLambda);
                }
                lambda[k] = newLambda;
            }
            updateR();
        }

        @Override
        public void updateAndQ(double k2, double temperature) {
            int n = lambda.length;
            double beta = 1.0/(temperature);
            double betaN = beta/n;
            double omegaN = Math.sqrt(n)/(beta*hbar);
            double omega2 = k2;
            for (int k=0; k<n; k++) {
                double s = Math.sin(Math.PI * k / n);
                lambda[k] = beta * 1 * (4.0 * omegaN * omegaN * s * s + omega2/n);
            }
            updateQ();
        }
    }

    public interface Transformer {

        /**
         * Uses the atom positions to update the given q
         */
        public void updateQ();

        /**
         * Updates internal data for the given temperature and k2, scales the given q values,
         * then updates real atom positions
         */
        public void update(double k2, double temperature);

        /**
         * Updates internal data for the given temperature and recomputes q for the current config
         */
        public void updateAndQ(double k2, double temperature);
    }


    private static class ModifierTemperature implements Modifier {
        private final int nBeads;
        private final SimQuantumAO sim;
        private final double hbar;

        public ModifierTemperature(int nBeads, SimQuantumAO sim, double hbar) {
            this.hbar = hbar;
            this.nBeads = nBeads;
            this.sim = sim;
        }

        @Override
        public void setValue(double newValue) {
            double newT = Math.exp(newValue);
            double beta = 1.0/(newT);
            double betaN = beta/ nBeads;
            double omegaN = Math.sqrt(nBeads)/(beta*hbar);

            sim.k2_kin = nBeads == 1 ? 0 : omegaN*omegaN;
            // p2Bond.setk2
            ((MCMoveHOReal) sim.move).setTemperature(newT);
            // doesn't really matter, but we just use this to hold T
            sim.integrator.setTemperature(newT);
        }

        @Override
        public double getValue() {
            return Math.log(sim.integrator.getTemperature());
        }

        @Override
        public Dimension getDimension() {
            return Null.DIMENSION;
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
                omega2 = omega2*(1.0 + omega2/12.0/ (sim.nBeads*omegaN*omegaN));
            }
            ((MCMoveHOReal) sim.move).setOmega2(omega2);

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

    private static class NewConfigAction implements IAction {

        private final SimQuantumAO sim;
        private final SimulationGraphic simGraphic;
        public Transformer transformer;

        public NewConfigAction(SimQuantumAO sim, SimulationGraphic simGraphic, Transformer transformer) {
            this.sim = sim;
            this.simGraphic = simGraphic;
            this.transformer = transformer;
        }

        public void setTransformer(Transformer transformer) {
            this.transformer = transformer;
        }

        @Override
        public void actionPerformed() {
            sim.move.doTrial();
            sim.move.acceptNotify();
            transformer.updateQ();
            simGraphic.getDisplayBox(sim.box).repaint();
        }
    }

    private static class FullUpdateAction implements IAction {
        private Transformer transformer;
        private final SimQuantumAO sim;
        private final SimulationGraphic simGraphic;

        public FullUpdateAction(SimQuantumAO sim, SimulationGraphic simGraphic, Transformer transformer) {
            this.transformer = transformer;
            this.sim = sim;
            this.simGraphic = simGraphic;
        }

        public void setTransformer(Transformer transformer) {
            this.transformer = transformer;
        }

        @Override
        public void actionPerformed() {
            transformer.update(sim.p1ah.getSpringConstants()[0], sim.integrator.getTemperature());
            simGraphic.getDisplayBox(sim.box).repaint();
        }
    }
}
