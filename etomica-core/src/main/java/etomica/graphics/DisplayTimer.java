/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountTime;
import etomica.integrator.IntegratorMD;

/**
 * DisplayBox to present the elapsed time in a
 * molecular dynamics simulation.  Designed for use with
 * a single integrator.
 */
public class DisplayTimer extends DisplayTextBox {

    public DisplayTimer(IntegratorMD integrator) {
        this(integrator, new DataSourceCountTime(integrator));
    }

    private DisplayTimer(IntegratorMD integrator, DataSourceCountTime timer) {
        super();
        putDataInfo(timer.getDataInfo());
        this.timer = timer;
        this.integrator = integrator;
        dataPump = new DataPumpListener(timer, this);
        integrator.getEventManager().addListener(dataPump);
        setUpdateInterval(10);
        setPrecision(7);
        graphic().setSize(100, 60);
    }

    /**
     * Sets the period for updating the display.  Number of integrator
     * interval events between updates.  Does not have any effect on the
     * value displayed; affects only how often it is updated.
     */
    public void setUpdateInterval(int interval) {
        dataPump.setInterval(interval);
    }

    /**
     * Unhooks the DisplayTimer from the integrator
     */
    public void dispose() {
        integrator.getEventManager().removeListener(dataPump);
    }

    /**
     * Returns the data source used to count the time, to permit
     * access to its methods for reset, etc.
     *
     * @return
     */
    public DataSourceCountTime getTimer() {
        return timer;
    }

    protected final DataPumpListener dataPump;
    protected final IntegratorMD integrator;
    protected final DataSourceCountTime timer;
}
