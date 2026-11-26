package etomica.modules.glass;

import etomica.atom.IAtom;
import etomica.data.ConfigurationStorage;
import etomica.data.meter.MeterStructureFactor;
import etomica.space.Vector;

/**
 * This allows mobility (or motion) structure factors to be computed taking the "position" of an atom to be the average
 * of the original and current positions of the atoms.
 */
public class AtomPositionConfig extends MeterStructureFactor.AtomPositionSource {
    protected final ConfigurationStorage configStorage;

    public AtomPositionConfig(ConfigurationStorage configStorage) {
        this.configStorage = configStorage;
    }

    public boolean ready() {
        return true;
    }

    public ConfigurationStorage getConfigStorage() {
        return configStorage;
    }

    public Vector position(IAtom atom) {
        int lastIndex = configStorage.getLastConfigIndex();
        if (lastIndex < 0) throw new RuntimeException("not enough configs to compute signal");
        Vector[] positions = configStorage.getSavedConfig(0);
        return positions[atom.getLeafIndex()];
    }
}
