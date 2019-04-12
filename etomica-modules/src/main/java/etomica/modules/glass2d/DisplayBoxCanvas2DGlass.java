package etomica.modules.glass2d;

import etomica.action.activity.Controller;
import etomica.atom.IAtom;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvas2D;
import etomica.space.Space;
import etomica.space.Vector;

import java.awt.*;

public class DisplayBoxCanvas2DGlass extends DisplayBoxCanvas2D {

    protected final ConfigurationStorage configStorage;
    protected int configIndex;
    protected final Vector dr;
    protected boolean drawDisplacement;

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
        double sigma = 0.5 * displayBox.getDiameterHash().getDiameter(a);
        // default diameter
        if (sigma == -1) sigma = 1;
        sigmaP = (int) (toPixels * sigma);
        sigmaP = (sigmaP == 0) ? 1 : sigmaP;
        xP = baseXP - (sigmaP >> 1);
        yP = baseYP - (sigmaP >> 1);
        g.fillOval(xP, yP, sigmaP, sigmaP);
        /* Draw the orientation line, if any */

        dr.Ev1Mv2(r, rOld);
        int dxy = (int) (toPixels * 0.5 * sigma);
        int dx = (int) (toPixels * dr.getX(0));
        int dy = (int) (toPixels * dr.getX(1));
        xP += dxy;
        yP += dxy;
        g.drawLine(xP, yP, xP - dx, yP - dy);
    }

}
