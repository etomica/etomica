package etomica.modules.glass;

import etomica.data.DataPump;

/**
 * Class that invokes the a pump only when an appropriate new configuration is available.
 */
public class ConfigurationStoragePumper implements ConfigurationStorage.ConfigurationStorageListener {
    private final ConfigurationStorage configStorage;
    private final DataPump pump;
    private int prevConfigIndex;

    public ConfigurationStoragePumper(DataPump pump, ConfigurationStorage configStorage) {
        this.configStorage = configStorage;
        this.pump = pump;
    }

    public void setPrevConfig(int index) {
        prevConfigIndex = index;
    }

    @Override
    public void newConfigruation() {
        long step = configStorage.getSavedSteps()[0];
        if (step==0 || step % (1L << prevConfigIndex) != 0) return;
        pump.actionPerformed();
    }
}
