package etomica.modules.glass;

import etomica.atom.IAtom;
import etomica.data.ConfigurationStorage;
import etomica.data.meter.MeterStructureFactor;
import etomica.space.Vector;

/**
 * This allows mobility (or motion) structure factors to be computed taking the "position" of an atom to be the average
 * of the original and current positions of the atoms.
 */
public class AtomPositionMobility extends MeterStructureFactor.AtomPositionSource {
    protected final ConfigurationStorage configStorage;
    protected final Vector dr;
    protected int prevConfigIndex;

    public AtomPositionMobility(ConfigurationStorage configStorage) {
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

    public ConfigurationStorage getConfigStorage() {
        return configStorage;
    }

    public Vector position(IAtom atom) {
        int idx = prevConfigIndex;
        int lastIndex = configStorage.getLastConfigIndex();
        if (lastIndex < idx) throw new RuntimeException("not enough configs to compute signal");
        Vector[] positions = configStorage.getSavedConfig(0);
        Vector[] prevPositions = configStorage.getSavedConfig(idx);
        int atomIndex = atom.getLeafIndex();
        dr.Ev1Pv2(positions[atomIndex], prevPositions[atomIndex]);
        dr.TE(0.5);
        return dr;
    }
}
