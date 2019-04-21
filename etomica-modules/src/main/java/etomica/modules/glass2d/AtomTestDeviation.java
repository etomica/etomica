package etomica.modules.glass2d;

import etomica.atom.AtomTest;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.space.Space;
import etomica.space.Vector;

public class AtomTestDeviation implements AtomTest {
    protected final ConfigurationStorage configStorage;
    protected int configIndex;
    protected final Box box;
    protected final Vector dr;
    protected double minDistance;

    public AtomTestDeviation(Box box, ConfigurationStorage configStorage) {
        this.box = box;
        this.configStorage = configStorage;
        Space space = box.getSpace();
        dr = space.makeVector();

        setMinDistance(0);
        configIndex = 100;
    }

    public void setConfigIndex(int idx) {
        configIndex = idx;
    }

    public double getMinDistance() { return this.minDistance;}

    public void setMinDistance(double minDistance) {
        this.minDistance = minDistance;
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

    @Override
    public boolean test(IAtom a) {
        return getDisplacementSq(a) > minDistance * minDistance;
    }
}
