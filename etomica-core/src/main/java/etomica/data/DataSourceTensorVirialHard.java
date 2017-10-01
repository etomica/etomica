/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;
import etomica.data.types.DataTensor;
import etomica.data.types.DataTensor.DataInfoTensor;
import etomica.integrator.IntegratorHard;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.units.dimensions.Energy;

/**
 * A MeterTensor that returns the virial component of the pressure tensor for a hard potential.  
 * This is the counterpart to MeterVelocityTensor, also called by the surface tension meter.
 *
 * @author Rob Riggleman
 */
public class DataSourceTensorVirialHard implements IDataSource, IntegratorHard.CollisionListener, java.io.Serializable {

    public DataSourceTensorVirialHard(Space space) {
        data = new DataTensor(space);
        dataInfo = new DataInfoTensor("Virial", Energy.DIMENSION, space);
        work = space.makeTensor();
        tag = new DataTag();
        dataInfo.addTag(tag);
    }

    public DataSourceTensorVirialHard(Space space, IntegratorHard integrator) {
        this(space);
        setIntegrator(integrator);
    }

    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    /**
     * Current value of the meter, obtained by dividing sum of collision virial contributions by time elapsed since last call.
     * If elapsed-time interval is zero, returns the value reported at the last call to the method.
     */
    public IData getData() {
        double currentTime = integratorHard.getCurrentTime();
        double elapsedTime = currentTime - lastTime;
        lastTime = currentTime;
        if(elapsedTime == 0.0) {
            data.E(Double.NaN);
            return data;
        }

        work.TE(-1./elapsedTime);
        data.x.E(work);
        work.E(0);
        return data;
    }

    /**
     * Sums contribution to virial for each collision.
     */
    public void collisionAction(IntegratorHard.Agent agent) {
        work.PE(agent.collisionPotential.lastCollisionVirialTensor());
    }

    /**
     * Contribution to the virial from the most recent collision of the given pair/potential.
     */
    public Tensor collisionValue(IntegratorHard.Agent agent) {
        data.x.E(agent.collisionPotential.lastCollisionVirialTensor());
        data.x.TE(1/(double)integratorHard.getBox().getLeafList().getAtomCount());
        return data.x;
    }

    /**
     * Informs meter of the integrator for its box, and zeros elapsed-time counter
     */
    public void setIntegrator(IntegratorHard newIntegrator) {
        if(newIntegrator == integratorHard) return;
        if(integratorHard != null) {
            integratorHard.removeCollisionListener(this);
        }
        integratorHard = newIntegrator;
        if(newIntegrator != null) {
            newIntegrator.addCollisionListener(this);
            lastTime = newIntegrator.getCurrentTime();
        }
    }

    public IntegratorHard getIntegrator() {
        return integratorHard;
    }
    
    private static final long serialVersionUID = 1L;
    protected double lastTime;
    protected IntegratorHard integratorHard;
    protected final DataTensor data;
    protected final Tensor work;
    protected final DataInfoTensor dataInfo;
    protected final DataTag tag;
}
