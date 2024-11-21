package etomica.modules.glass;

import etomica.data.ConfigurationStorage;
import etomica.data.DataPump;

/**
 * Class that invokes the a pump only when an appropriate new configuration is available.
 *
 * With only prevStep set, the pump is fired every 2^prevStep steps; if prevStep=4, then
 * the pump is fired for step 0, 16, 32, 48, 64...
 *
 * With bigStep also set (and bigStep>prevStep), then the pump will fire every 2^bigSteps
 * steps; if bigStep=4, then the pump is fired for step 0, 16, 32, 48, 64.... It is also
 * fired prevStep steps before that; if bigStep=4 and prevStep=1, then the pump will fire
 * for step 0, 14, 16, 30, 32, 46, 48, 62, 64...
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
        if (step>0 && ((step % (1L << Math.max(bigStep, prevStep)) == 0) ||
                (bigStep > prevStep && ((step + (1L << prevStep)) % (1L << bigStep) == 0)))) {
            pump.actionPerformed();
        }
    }
}
