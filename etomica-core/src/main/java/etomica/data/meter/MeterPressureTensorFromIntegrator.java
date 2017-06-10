/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.data.DataInfo;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataTensor;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.space.Space;
import etomica.units.dimensions.Pressure;

/**
 * Acts as a DataSource to retrieve the pressure from the integrator.
 * Currently only works with IntegratorVelocityVerlet (but could work just as
 * well with IntegratorVerlet;
 * see https://rheneas.eng.buffalo.edu/bugzilla/show_bug.cgi?id=164 )
 */
public class MeterPressureTensorFromIntegrator implements IEtomicaDataSource, java.io.Serializable {

    public MeterPressureTensorFromIntegrator(Space space) {
        tag = new DataTag();
        this.space = space;
    }
    
    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    public void setIntegrator(IntegratorVelocityVerlet newIntegrator) {
        integrator = newIntegrator;
        if (integrator != null) {
            data = new DataTensor(space);
            dataInfo = new DataTensor.DataInfoTensor("Pressure", Pressure.DIMENSION, space);
        }
    }
    
    public IntegratorBox getIntegrator() {
        return integrator;
    }
    
    public IData getData() {
        data.x.E(integrator.getPressureTensor());
        return data;
    }
    
    private static final long serialVersionUID = 1L;
    protected IntegratorVelocityVerlet integrator;
    protected DataTensor data;
    protected DataInfo dataInfo;
    protected final DataTag tag;
    private final Space space;
}
