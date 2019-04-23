/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.data.DataSourceScalar;
import etomica.units.dimensions.Time;

public class DataSourcePrevTime extends DataSourceScalar {
    protected ConfigurationStorage configStorage;
    protected int prevConfigIndex;

    public DataSourcePrevTime(ConfigurationStorage configStorage) {
        super("previous time", Time.DIMENSION);
        this.configStorage = configStorage;
    }

    public void setPrevConfigIndex(int idx) {
        prevConfigIndex = idx;
    }

    public void setConfigStorage(ConfigurationStorage configStorage){
        this.configStorage = configStorage;
    }

    public ConfigurationStorage getConfigStorage(){
        return configStorage;
    }

    @Override
    public double getDataAsScalar() {
        int lastIndex = configStorage.getLastConfigIndex();
        if (lastIndex < 1) return 0;
        int myIndex = prevConfigIndex;
        if (myIndex > lastIndex) myIndex = lastIndex;
        double prevTime = configStorage.getSavedTimes()[myIndex];
        double nowTime = configStorage.getSavedTimes()[0];
        return nowTime - prevTime;
    }
}
