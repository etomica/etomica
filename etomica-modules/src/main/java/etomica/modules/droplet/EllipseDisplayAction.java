/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.droplet;

import java.awt.Color;

import etomica.simulation.Simulation;
import org.jmol.util.Point3f;

import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.SimulationGraphic;
import g3dsys.control.G3DSys;
import g3dsys.images.Ellipse;

public class EllipseDisplayAction {
    public EllipseDisplayAction(SimulationGraphic graphic, float dropDiameter) {
        this.graphic = graphic;
        nominalDiameter = dropDiameter;
        Simulation sim = graphic.getSimulation();
        G3DSys gsys = ((DisplayBoxCanvasG3DSys)graphic.getDisplayBox(sim.getBox(0)).canvas).getG3DSys();
        ellipse = new Ellipse(gsys, G3DSys.getColix(Color.WHITE), Point3f.new3(0,0,3), nominalDiameter, 1);
    }
    public synchronized void displayEllipse(double g) {
        final Simulation sim = graphic.getSimulation();
        if (g == oldG) {
            return;
        }
        oldG = g;
        final G3DSys gsys = ((DisplayBoxCanvasG3DSys)graphic.getDisplayBox(sim.getBox(0)).canvas).getG3DSys();
        gsys.rotateToHome();
        gsys.rotateByY(90);
  
        double d = def[(int)Math.round(g*10)];
        double factor = (1+d) / (1-d);
        float xdiam = (float)(Math.pow(factor, -1.0/3.0));
        float zdiam = 1/(xdiam*xdiam);
        float aspect = xdiam/zdiam;
        xdiam *= nominalDiameter;
    
        ellipse.setD(2.2f*xdiam);
        ellipse.setAspect(aspect);
        serial++;
        if (!displayed) {
            gsys.addFig(ellipse);
            displayed = true;
        }

        if (undisplayThread != null) {
            synchronized (undisplayThread) {
                undisplayThread.notify();
            }
        }
        undisplayThread = new Thread(){
            public void run() {
                int oldSerial = serial;
                try {
                    synchronized(this) {
                        wait(1000);
                    }
                }
                catch (InterruptedException e) {}
                synchronized (this) {
                    if (serial == oldSerial) {
                        // we were the last one to mess with the ellipse, so we can remove it
                        gsys.removeFig(ellipse);
                        graphic.getDisplayBox(sim.getBox(0)).repaint();
                        displayed = false;
                    }
                    if (this == undisplayThread) {
                        //dispose of ourselves!
                        undisplayThread = null;
                    }
                }
            }
        };
        undisplayThread.start();
        ellipse.setX((float)(-sim.getBox(0).getBoundary().getBoxSize().getX(0)*0.6));
        ellipse.setY(0);
        ellipse.setZ(0);
        graphic.getDisplayBox(sim.getBox(0)).repaint();
    }

    protected Thread undisplayThread;
    protected final SimulationGraphic graphic;
    protected final Ellipse ellipse;
    protected double oldG;
    protected boolean displayed = false;
    protected int serial;
    protected float nominalDiameter;
    protected double[] def = {0.000000, 0.011460, 0.022439, 0.032981, 0.043126, 0.052907, 0.062354,
                              0.071492, 0.080344, 0.088929, 0.097265, 0.105368, 0.113251, 0.120929,
                              0.128412, 0.135711, 0.142836, 0.149795, 0.156596, 0.163248, 0.169756,
                              0.176127, 0.182368, 0.188482, 0.194477, 0.200356, 0.206124, 0.211784,
                              0.217342, 0.222800, 0.228162, 0.233431, 0.238611, 0.243704, 0.248713,
                              0.253641, 0.258489, 0.263261, 0.267959, 0.272585, 0.277140, 0.281627,
                              0.286048, 0.290404, 0.294697, 0.298929, 0.303100, 0.307214, 0.311271,
                              0.315272, 0.319219, 0.323113, 0.326955, 0.330747, 0.334489, 0.338183,
                              0.341830, 0.345430, 0.348985, 0.352496, 0.355963, 0.359388, 0.362771,
                              0.366113, 0.369415, 0.372678, 0.375902, 0.379088, 0.382238, 0.385351,
                              0.388428, 0.391470, 0.394478, 0.397452, 0.400393, 0.403302, 0.406178,
                              0.409023, 0.411837, 0.414620, 0.417374, 0.420098, 0.422794, 0.425461,
                              0.428100, 0.430711, 0.433296, 0.435854, 0.438385, 0.440891, 0.443372,
                              0.445828, 0.448259, 0.450666, 0.453048, 0.455408, 0.457745, 0.460058,
                              0.462350, 0.464619, 0.466866};
}

