package etomica.modules.glass;

import etomica.atom.IAtom;
import etomica.data.meter.MeterStructureFactor;
import etomica.space.Vector;

public class AtomSignalMotion extends MeterStructureFactor.AtomSignalSourceByType {
    protected final ConfigurationStorage configStorage;
    protected int prevConfigIndex;
    protected int xyz;

    public AtomSignalMotion(ConfigurationStorage configStorage, int xyz) {
        this.configStorage = configStorage;
        this.xyz = xyz;
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
        int idx = prevConfigIndex;
        int lastIndex = configStorage.getLastConfigIndex();
        if (lastIndex < idx) throw new RuntimeException("not enough configs to compute signal");
        Vector[] positions = configStorage.getSavedConfig(0);
        Vector[] prevPositions = configStorage.getSavedConfig(idx);
        int atomIndex = atom.getLeafIndex();
        double dxyz = positions[atomIndex].getX(xyz) - prevPositions[atomIndex].getX(xyz);
        return s * dxyz;
    }
}
