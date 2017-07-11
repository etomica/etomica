/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.rosmosis;

import etomica.box.Box;
import etomica.space.Vector;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.units.dimensions.Pressure;

/**
 * Meter to measure flux across a boundary or boundaries.  If an atom is on one
 * side of the boundary at one time and on the other side of the boundary the
 * next time, then the meter counts that as a crossing.  Each boundary can have
 * a coefficient, such that flow into one region through boundaries on either
 * side can be considered "positive" flux.
 * 
 * If an IntegratorMD is used, flux will be given in terms of crossings per
 * area per time.  Otherwise, flux will be given in terms of crossings per
 * area per step.
 *
 * @author Andrew Schultz
 */
public class MeterOsmoticPressure implements IEtomicaDataSource {

    public MeterOsmoticPressure(IPotentialCalculationWallForce pc, Box box) {
        this.pc = pc;
        this.box = box;
        data = new DataDouble();
        dataInfo = new DataInfoDouble("osmotic pressure", Pressure.DIMENSION);
        tag = new DataTag();
    }
    
    public IData getData() {
        Vector dimensions = box.getBoundary().getBoxSize();
        data.x = -pc.getWallForce() / (dimensions.getX(1) * dimensions.getX(2));
        return data;
    }

    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    protected final IPotentialCalculationWallForce pc;
    protected final Box box;
    protected DataDouble data;
    protected DataInfoDouble dataInfo;
    protected DataTag tag;
}
