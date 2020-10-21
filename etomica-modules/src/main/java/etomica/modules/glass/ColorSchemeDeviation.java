/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.graphics.ColorScheme;
import etomica.space.Space;
import etomica.space.Vector;

import java.awt.*;

public class ColorSchemeDeviation extends ColorScheme {
    protected ConfigurationStorage configStorage;
    protected int configIndex;
    protected final Box box;
    protected final Vector dr;
    protected final Color[] colors;
    protected double fac;
    protected boolean isString;
    protected double dr2StringMax = 0.6*0.6;


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

    public void setIsString(boolean isString) {this.isString = isString;}

    public void setDrString(double drString) {this.dr2StringMax = drString*drString;}

    public double getDrString() {return Math.sqrt(this.dr2StringMax);}

    public boolean getIsString(){return isString;}

    public int getConfigIndex() {
        return configIndex;
    }

    public void setConfigIndex(int idx) {
        configIndex = idx;
    }

    public void setLengthFactor(double lenghtFactor) {
        this.fac = lenghtFactor;
    }

    public void setConfigStorage(ConfigurationStorage configStorage){
        this.configStorage = configStorage;
    }

    public ConfigurationStorage getConfigStorage(){return configStorage;}

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

    public boolean isAtomInString(IAtom atom){
        int idx = configIndex;
        int lastIndex = configStorage.getLastConfigIndex();
        if (idx > lastIndex) idx = lastIndex;
        if (idx == -1) return false;
        Vector r = configStorage.getSavedConfig(0)[atom.getLeafIndex()];
        for (int j = 0; j < box.getLeafList().size(); j++) {
            if (atom.getLeafIndex() == j) continue;
            Vector oldR = configStorage.getSavedConfig(idx)[j];
            dr.Ev1Mv2(r,oldR);
            box.getBoundary().nearestImage(dr);
            double dr2 = dr.squared();
            if(dr2 < dr2StringMax){
                return true;
            }
        }
        return false;
    }

    @Override
    public Color getAtomColor(IAtom a) {
        if(isString && isAtomInString(a)){
            return dr.getD()==3?Color.WHITE:Color.BLACK;
        }
        double s = getRelativeDisplacement(a);
        if (s >= 2) return colors[510];
        int i = (int) (510.9999 * s / 2);
        return colors[i];
    }
}
