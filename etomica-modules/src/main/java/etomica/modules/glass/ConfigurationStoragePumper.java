package etomica.modules.glass;

import etomica.data.DataPump;

/**
 * Class that invokes the a pump only when an appropriate new configuration is available.
 */
public class ConfigurationStoragePumper implements ConfigurationStorage.ConfigurationStorageListener {
    private final ConfigurationStorage configStorage;
    private final DataPump pump;
    private int prevStep;
    private int bigStep;

    public ConfigurationStoragePumper(DataPump pump, ConfigurationStorage configStorage) {
        this.configStorage = configStorage;
        this.pump = pump;
    }

    public void setBigStep(int bigStep) {
        this.bigStep = bigStep;
    }

    public void setPrevStep(int step) {
        prevStep = step;
    }

    @Override
    public void newConfigruation() {
        long step = configStorage.getSavedSteps()[0];
        if ((step % (1L << Math.max(bigStep, prevStep)) == 0) ||
                (bigStep > prevStep && ((step + (1L << prevStep)) % (1L << bigStep) == 0))) {
            pump.actionPerformed();
        }
    }
}
