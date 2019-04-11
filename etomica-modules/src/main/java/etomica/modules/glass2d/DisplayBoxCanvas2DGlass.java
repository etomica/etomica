package etomica.modules.glass2d;

import etomica.action.activity.Controller;
import etomica.atom.IAtom;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvas2D;
import etomica.space.Space;
import etomica.space.Vector;

import java.awt.*;
import java.util.ArrayList;

public class DisplayBoxCanvas2DGlass extends DisplayBoxCanvas2D {

    protected final java.util.List<Vector> oldPositions = new ArrayList<>();
    protected final Vector dr;

    public DisplayBoxCanvas2DGlass(DisplayBox _box, Space _space, Controller controller) {
        super(_box, _space, controller);
        reset();
        dr = _space.makeVector();
    }

    public void reset() {
        oldPositions.clear();
        for (IAtom a : displayBox.getBox().getLeafList()) {
            Vector v = space.makeVector();
            v.E(a.getPosition());
            oldPositions.add(v);
        }
    }

    protected void drawAtom(Graphics g, int[] origin, IAtom a) {
        Vector r = a.getPosition();
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

        dr.Ev1Mv2(r, oldPositions.get(a.getLeafIndex()));
        displayBox.getBox().getBoundary().nearestImage(dr);
        int dxy = (int) (toPixels * 0.5 * sigma);
        int dx = (int) (toPixels * dr.getX(0));
        int dy = (int) (toPixels * dr.getX(1));
        xP += dxy;
        yP += dxy;
        g.drawLine(xP, yP, xP - dx, yP - dy);
    }

}
