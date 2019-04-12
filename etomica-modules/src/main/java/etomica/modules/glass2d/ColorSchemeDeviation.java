package etomica.modules.glass2d;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.graphics.ColorScheme;
import etomica.space.Space;
import etomica.space.Vector;

import java.awt.*;

public class ColorSchemeDeviation extends ColorScheme {
    protected final ConfigurationStorage configStorage;
    protected int configIndex;
    protected final Box box;
    protected final Vector dr;
    protected final Color[] colors;
    protected double fac;

    public ColorSchemeDeviation(Box box, ConfigurationStorage configStorage) {
        this.box = box;
        this.configStorage = configStorage;
        Space space = box.getSpace();
        dr = space.makeVector();
        colors = new Color[511];
        for (int i = 0; i < 256; i++) {
            colors[i] = new Color(0, i, 255 - i);
        }
        for (int i = 1; i < 256; i++) {
            colors[255 + i] = new Color(i, 255 - i, 0);
        }
        setLengthFactor(1);
        configIndex = 100;
    }

    public void setConfigIndex(int idx) {
        configIndex = idx;
    }

    public void setLengthFactor(double lenghtFactor) {
        this.fac = lenghtFactor;
    }

    public double getDisplacementSq(IAtom a) {
        int idx = configIndex;
        int lastIndex = configStorage.getLastConfigIndex();
        if (idx > lastIndex) idx = lastIndex;
        if (idx == -1) return 0;
        Vector r = configStorage.getSavedConfig(0)[a.getLeafIndex()];
        Vector oldR = configStorage.getSavedConfig(idx)[a.getLeafIndex()];

        dr.Ev1Mv2(r, oldR);
        return dr.squared();
    }

    public double getRelativeDisplacement(IAtom a) {
        return Math.sqrt(getDisplacementSq(a)) / fac;
    }

    @Override
    public Color getAtomColor(IAtom a) {
        double s = getRelativeDisplacement(a);
        if (s >= 2) return colors[510];
        int i = (int) (510.9999 * s / 2);
        return colors[i];
    }
}
