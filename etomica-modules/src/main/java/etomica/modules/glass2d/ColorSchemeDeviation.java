package etomica.modules.glass2d;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.graphics.ColorScheme;
import etomica.space.Space;
import etomica.space.Vector;

import java.awt.*;

public class ColorSchemeDeviation extends ColorScheme {
    protected final Vector[] scaledCoords0;
    protected final Box box;
    protected final Vector dr;
    protected final Color[] colors;
    protected double fac;

    public ColorSchemeDeviation(Box box) {
        this.box = box;
        Space space = box.getSpace();
        dr = space.makeVector();
        colors = new Color[511];
        for (int i = 0; i < 256; i++) {
            colors[i] = new Color(0, i, 255 - i);
        }
        for (int i = 1; i < 256; i++) {
            colors[255 + i] = new Color(i, 255 - i, 0);
        }
        int numAtoms = box.getLeafList().size();
        scaledCoords0 = space.makeVectorArray(numAtoms);
        setLengthFactor(1);
        reset();
    }

    public void reset() {
        IAtomList atoms = box.getLeafList();
        for (IAtom a : atoms) {
            scaledCoords0[a.getLeafIndex()].E(a.getPosition());
        }
    }

    public void setLengthFactor(double lenghtFactor) {
        this.fac = lenghtFactor;
    }

    public Vector getDisplacement(IAtom a) {
        dr.Ev1Mv2(a.getPosition(), scaledCoords0[a.getLeafIndex()]);
        box.getBoundary().nearestImage(dr);
        return dr;
    }

    public double getRelativeDisplacement(IAtom a) {
        return Math.sqrt(getDisplacement(a).squared()) / fac;
    }

    @Override
    public Color getAtomColor(IAtom a) {
        double s = getRelativeDisplacement(a);
        if (s >= 2) return colors[510];
        int i = (int) (510.9999 * s / 2);
        return colors[i];
    }
}
