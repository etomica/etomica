/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation;

import etomica.data.IEtomicaDataSource;

public class DataStreamHeader implements java.io.Serializable {

    private final IEtomicaDataSource dataSource;
    private Object client;

    DataStreamHeader(IEtomicaDataSource dataSource, Object object) {
        this.dataSource = dataSource;
        client = object;
    }

    public IEtomicaDataSource getDataSource() {
        return dataSource;
    }

    public Object getClient() {
        return client;
    }

    public void setClient(Object object) {
        client = object;
    }
}
