package etomica.modules.glass;

import etomica.atom.IAtom;
import etomica.data.ConfigurationStorage;
import etomica.data.meter.MeterStructureFactor;
import etomica.space.Vector;

public class AtomSignalMobility extends MeterStructureFactor.AtomSignalSourceByType {
    protected final ConfigurationStorage configStorage;
    protected final Vector dr;
    protected int prevConfigIndex;

    public AtomSignalMobility(ConfigurationStorage configStorage) {
        this.configStorage = configStorage;
        dr = configStorage.getBox().getSpace().makeVector();
    }

    public boolean ready() {
        int lastIndex = configStorage.getLastConfigIndex();
        return prevConfigIndex <= lastIndex;
    }

    public void setPrevConfig(int prevConfigIndex) {
        this.prevConfigIndex = prevConfigIndex;
    }

    public int getPrevConfigIndex() {
        return prevConfigIndex;
    }

    public double signal(IAtom atom) {
        double s = super.signal(atom);
        if (s == 0) return 0;
        int idx = prevConfigIndex;
        int lastIndex = configStorage.getLastConfigIndex();
        if (lastIndex < idx) throw new RuntimeException("not enough configs to compute signal");
        Vector[] positions = configStorage.getSavedConfig(0);
        Vector[] prevPositions = configStorage.getSavedConfig(idx);
        int atomIndex = atom.getLeafIndex();
        dr.Ev1Mv2(positions[atomIndex], prevPositions[atomIndex]);
        return s * Math.sqrt(dr.squared());
    }
}
