package etomica.modules.glass2d;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.graphics.ColorScheme;
import etomica.space.Space;
import etomica.space.Vector;

import java.awt.*;

public class ColorSchemeDirection extends ColorScheme {
    protected ConfigurationStorage configStorage;
    protected int configIndex;
    protected final Box box;
    protected final Vector dr;
    protected final Color[] colors;
    protected double fac;
    protected int axis;

    public ColorSchemeDirection(Box box, ConfigurationStorage configStorage) {
        this.box = box;
        this.configStorage = configStorage;
        Space space = box.getSpace();
        dr = space.makeVector();
        colors = new Color[255 * 6 + 1];
        int j = 0;
        for (int i = 0; i < 256; i++) {
            colors[j++] = new Color(0, i, 255);
        }
        for (int i = 254; i >= 0; i--) {
            colors[j++] = new Color(0, 255, i);
        }
        for (int i = 1; i < 256; i++) {
            colors[j++] = new Color(i, 255, 0);
        }
        for (int i = 254; i >= 0; i--) {
            colors[j++] = new Color(255, i, 0);
        }
        for (int i = 1; i < 256; i++) {
            colors[j++] = new Color(255, 0, i);
        }
        for (int i = 254; i >= 0; i--) {
            colors[j++] = new Color(i, 0, 255);
        }
        configIndex = 100;
        setAxis(2);
    }

    public void setAxis(int i) {
        axis = i;
    }

    public int getAxis() {
        return axis;
    }

    public int getConfigIndex() {
        return configIndex;
    }

    public void setConfigIndex(int idx) {
        configIndex = idx;
    }

    public void setConfigStorage(ConfigurationStorage configStorage){
        this.configStorage = configStorage;
    }

    public ConfigurationStorage getConfigStorage(){ return configStorage;}

    public Vector getDisplacement(IAtom a) {
        int idx = configIndex;
        int lastIndex = configStorage.getLastConfigIndex();
        if (idx > lastIndex) idx = lastIndex;
        if (idx == -1) {
            dr.E(0);
            return dr;
        }
        Vector r = configStorage.getSavedConfig(0)[a.getLeafIndex()];
        Vector oldR = configStorage.getSavedConfig(idx)[a.getLeafIndex()];

        dr.Ev1Mv2(r, oldR);
        return dr;
    }

    @Override
    public Color getAtomColor(IAtom a) {
        double[] x = new double[2];
        int j = 0;
        getDisplacement(a);
        for (int i = 0; i < dr.getD(); i++) {
            if (i == axis) continue;
            x[j] = dr.getX(i);
            j++;
        }
        if (x[0] == 0.0 && x[1] == 0.0) return Color.BLACK;
        double theta = Math.atan2(x[1], x[0]);
        if (theta < 0) theta += 2 * Math.PI;
        double s = theta / (2 * Math.PI);
        if (s >= 1) s = 1;
        int i = (int) ((colors.length - 0.00001) * s);
        return colors[i];
    }
}
