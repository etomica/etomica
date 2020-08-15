/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode.nptdemo;
import javax.swing.JFrame;

import etomica.action.IAction;
import etomica.atom.DiameterHashByType;
import etomica.graphics.ColorScheme;
import etomica.graphics.ColorSchemeOverlap;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DeviceToggleButton;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayCanvas;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.modifier.ModifierBoolean;

/**
 * Simple hard-sphere molecular dynamics simulation in 2D.
 *
 * @author David Kofke
 */
 
public class HSNPT2DGraphic extends SimulationGraphic {
    
    private static final long serialVersionUID = 1L;
    public double nominalScale;
    public DisplayBox displayBoxScaled;
    public HSNPT2DSim sim;
    protected final ColorScheme redColorScheme;
    protected final ColorSchemeOverlap colorSchemeOverlap;
    protected final ColorSchemeScaledOverlap colorSchemeScaled;

    public HSNPT2DGraphic(HSNPT2DSim sim) {
        super(sim, TABBED_PANE);
        this.sim = sim;
        getDisplayBox(sim.box).setScale(0.8);
        redColorScheme = getDisplayBox(sim.box).getColorScheme();

        displayBoxScaled = new DisplayBox(sim, sim.box);
        displayBoxScaled.setScale(0.9);
        System.out.println("displayBoxScaled created");
        DisplayCanvas oldCanvas = displayBoxScaled.canvas;
        displayBoxScaled.setBoxCanvas(new DisplayBoxCanvas2DNpTScaling(displayBoxScaled, sim.getSpace(), sim.getController(), sim.coordinateDefinition));
        oldCanvas.dispose();
        System.out.println("displayBoxScaled has a new canvas");
        addAsTab(displayBoxScaled.graphic(), "scaled", true);
        IntegratorListenerAction repaintAction = new IntegratorListenerAction(new IAction() {
            public void actionPerformed() {
                displayBoxScaled.repaint();
            }
        });
        repaintAction.setInterval(100);
        sim.integrator.getEventManager().addListener(repaintAction);
        displayBoxScaled.setDiameterHash(getDisplayBox(sim.box).getDiameterHash());
        
        colorSchemeOverlap = new ColorSchemeOverlap(sim.getSpace(), sim.potentialMaster, sim.box);
        getDisplayBox(sim.box).setColorScheme(colorSchemeOverlap);
        
        colorSchemeScaled = new ColorSchemeScaledOverlap(sim.getSpace(), sim.potentialMaster, sim.coordinateDefinition);
        displayBoxScaled.setColorScheme(colorSchemeScaled);

        getController().getReinitButton().setPostAction(new IAction() {
            public void actionPerformed() {
                getDisplayBox(simulation.getBox(0)).repaint();
                displayBoxScaled.repaint();
            }
        });
        
        nominalScale = getDisplayBox(sim.box).getScale();

        DeviceSlider rhoSlider = new DeviceSlider(sim.getController(), this, "displayDensity");
        rhoSlider.setPrecision(3);
        rhoSlider.setMinimum(0.85);
        rhoSlider.setMaximum(1.15);
        rhoSlider.setNMajor(2);
        rhoSlider.setNMinor(3);
        rhoSlider.doUpdate();
        rhoSlider.setShowBorder(true);
        rhoSlider.setShowValues(true);
        rhoSlider.setLabel("density");
        add(rhoSlider);
        
        setPressure(15);
        
        DeviceSlider pSlider = new DeviceSlider(sim.getController(), this, "pressure");
        pSlider.setPrecision(1);
        pSlider.setMinimum(10);
        pSlider.setMaximum(20);
        pSlider.setNMajor(4);
        pSlider.doUpdate();
        pSlider.setShowBorder(true);
        pSlider.setShowValues(true);
        pSlider.setLabel("pressure");
        add(pSlider);

        DeviceToggleButton colorButton = new DeviceToggleButton(sim.getController(), new ModifierBoolean() {
            
            public void setBoolean(boolean b) {
                if (b) {
                    getDisplayBox(simulation.getBox(0)).setColorScheme(redColorScheme);
                    displayBoxScaled.setColorScheme(redColorScheme);
                }
                else {
                    getDisplayBox(simulation.getBox(0)).setColorScheme(colorSchemeOverlap);
                    displayBoxScaled.setColorScheme(colorSchemeScaled);
                }
            }
            
            public boolean getBoolean() {
                return displayBoxScaled.getColorScheme() == redColorScheme;
            }
        }, "Color: overlaps", "Color: red");
        add(colorButton);
    }
    
    public void setDisplayDensity(double newDisplayDensity) {
        double density = sim.box.getMoleculeList().size() / sim.box.getBoundary().volume();
        double newScale = Math.sqrt(density/newDisplayDensity) * nominalScale;
        DisplayBox unscaledDisplayBox = getDisplayBox(sim.box);
        unscaledDisplayBox.setScale(newScale);
        displayBoxScaled.setScale(newScale);
        ((DisplayBoxCanvas2DNpTScaling)displayBoxScaled.canvas).setDisplayDensity(newDisplayDensity);
        colorSchemeScaled.setDisplayDensity(newDisplayDensity);
        double displayDiameter = nominalScale / newScale;
        ((DiameterHashByType)unscaledDisplayBox.getDiameterHash()).setDiameter(sim.getSpecies(0).getAtomType(0), displayDiameter);
        colorSchemeOverlap.setSigma(displayDiameter);
        unscaledDisplayBox.repaint();
        displayBoxScaled.repaint();
    }
    
    public double getDisplayDensity() {
        double scale = getDisplayBox(sim.box).getScale();
        double density = sim.box.getMoleculeList().size() / sim.box.getBoundary().volume();
        double s = scale/nominalScale;
        return density /(s*s); 
    }
    
    public void setPressure(double newPressure) {
        sim.pressure = newPressure;
        ((DisplayBoxCanvas2DNpTScaling)displayBoxScaled.canvas).setPressure(newPressure);
        colorSchemeScaled.setPressure(newPressure);
        displayBoxScaled.repaint();
    }
    
    public double getPressure() {
        return sim.pressure;
    }

    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        HSNPT2DSim sim = new HSNPT2DSim();
        sim.getController().setSleepPeriod(1);
        HSNPT2DGraphic graphic = new HSNPT2DGraphic(sim);
        JFrame f = graphic.makeAndDisplayFrame();
        f.setSize(700, 500);

//        MeterPressureHard meterP = new MeterPressureHard(sim.space);
//        meterP.setIntegrator(sim.integrator);
//        sim.integrator.addCollisionListener(meterP);
//        sim.ai.setMaxSteps(10000);
//        sim.getController().actionPerformed();
//        System.out.println("P: "+meterP.getDataAsScalar());
//        System.exit(1);

   }
}
