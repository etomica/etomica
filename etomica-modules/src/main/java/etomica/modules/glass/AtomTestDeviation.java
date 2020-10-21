/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.atom.AtomTest;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.space.Space;
import etomica.space.Vector;

public class AtomTestDeviation implements AtomTest {
    protected ConfigurationStorage configStorage;
    protected int configIndex;
    protected final Box box;
    protected final Vector dr;
    protected double minDistance;
    protected boolean doMobileOnly = true;

    public AtomTestDeviation(Box box, ConfigurationStorage configStorage) {
        this.box = box;
        this.configStorage = configStorage;
        Space space = box.getSpace();
        dr = space.makeVector();

        setMinDistance(0);
        configIndex = 100;
    }

    public void setDoMobileOnly(boolean doMobileOnly) {
        this.doMobileOnly = doMobileOnly;
    }

    public boolean getDoMobileOnly() {
        return doMobileOnly;
    }

    public void setConfigIndex(int idx) {
        configIndex = idx;
    }

    public double getMinDistance() { return this.minDistance;}

    public void setMinDistance(double minDistance) {
        this.minDistance = minDistance;
    }

    public void setConfigStorage(ConfigurationStorage configStorage){
        this.configStorage = configStorage;
    }

    public ConfigurationStorage getConfigStorage(){ return configStorage; }

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

    @Override
    public boolean test(IAtom a) {
        return getDisplacementSq(a) > minDistance * minDistance == doMobileOnly;
    }
}
