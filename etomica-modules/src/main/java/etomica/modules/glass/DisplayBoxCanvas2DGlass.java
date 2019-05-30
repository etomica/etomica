/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.action.activity.Controller;
import etomica.atom.IAtom;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvas2D;
import etomica.space.Space;
import etomica.space.Vector;

import java.awt.*;

public class DisplayBoxCanvas2DGlass extends DisplayBoxCanvas2D implements DisplayBoxCanvasGlass {

    protected ConfigurationStorage configStorage;
    protected int configIndex;
    protected final Vector dr;
    protected boolean drawDisplacement, flipDisplacement;

    public DisplayBoxCanvas2DGlass(DisplayBox _box, Space _space, Controller controller, ConfigurationStorage configStorage) {
        super(_box, _space, controller);
        dr = _space.makeVector();
        this.configStorage = configStorage;
        configIndex = 100;
    }

    public void setConfigIndex(int idx) {
        configIndex = idx;
        repaint();
    }

    public int getConfigIndex() {
        return configIndex;
    }

    public void setDrawDisplacement(boolean doDrawDisplacement) {
        this.drawDisplacement = doDrawDisplacement;
    }

    public boolean getDrawDisplacement() {
        return drawDisplacement;
    }

    public void setConfigStorage(ConfigurationStorage configStorage){
        this.configStorage = configStorage;
    }

    public ConfigurationStorage getConfigStorage(){ return configStorage;}


    protected void drawAtom(Graphics g, int[] origin, IAtom a) {
        if (!drawDisplacement) {
            super.drawAtom(g, origin, a);
            return;
        }
        int idx = configIndex;
        int lastIndex = configStorage.getLastConfigIndex();
        if (idx > lastIndex) idx = lastIndex;
        Vector[] oldPositions = idx == -1 ? null : configStorage.getSavedConfig(idx);
        Vector r = idx == -1 ? a.getPosition() : configStorage.getSavedConfig(0)[a.getLeafIndex()];
        Vector rOld = idx == -1 ? r : configStorage.getSavedConfig(idx)[a.getLeafIndex()];
        int sigmaP, xP, yP, baseXP, baseYP;

        g.setColor(displayBox.getColorScheme().getAtomColor(a));

        double toPixels = pixel.toPixels() * displayBox.getScale();

        baseXP = origin[0] + (int) (toPixels * r.getX(0));
        baseYP = origin[1] + (int) (toPixels * r.getX(1));
        /* Draw the core of the atom, specific to the dimension */
        double sigma = displayBox.getDiameterHash().getDiameter(a);
        // default diameter
        if (sigma == -1) sigma = 1;
        sigmaP = (int) (toPixels * sigma);
        sigmaP = (sigmaP == 0) ? 1 : sigmaP;
        xP = baseXP - (sigmaP >> 1);
        yP = baseYP - (sigmaP >> 1);

        dr.Ev1Mv2(r, rOld);
        int dxy = (int) (toPixels * 0.5 * sigma);
        int dx = (int) (toPixels * dr.getX(0));
        int dy = (int) (toPixels * dr.getX(1));
        int xP1 = xP + dxy;
        int yP1 = yP + dxy;
        g.drawLine(xP1, yP1, xP1 - dx, yP1 - dy);

        if (flipDisplacement) {
            xP -= dx;
            yP -= dy;
        }
        g.fillOval(xP, yP, sigmaP, sigmaP);

    }

    public void setFlipDisplacement(boolean flipDisplacement) {
        this.flipDisplacement = flipDisplacement;
    }

    public boolean getFlipDisplacement() {
        return flipDisplacement;
    }

}
