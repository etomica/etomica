/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataTensor;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationPressureTensor;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.units.dimensions.Pressure;

/**
 * Meter for evaluation of the soft-potential pressure tensor in a box.  This
 * should only be used when using an MC Integrator.  For MD, use a
 * MeterPressureTensorFromIntegrator.
 *
 * @author Andrew Schultz
 */
public class MeterPressureTensor implements IEtomicaDataSource {
    
    public MeterPressureTensor(PotentialMaster potentialMaster, Space space) {
    	super();
        this.potentialMaster = potentialMaster;
        data = new DataTensor(space);
        tag = new DataTag();
        dataInfo = new DataTensor.DataInfoTensor("Pressure",Pressure.dimension(space.D()), space);
        dataInfo.addTag(tag);
        rD = 1.0/space.D();
        iteratorDirective = new IteratorDirective();
        iteratorDirective.includeLrc = true;
        pc = new PotentialCalculationPressureTensor(space);
    }

    public DataTag getTag() {
        return tag;
    }
    
    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }
    
    /**
     * Sets the integrator associated with this instance.  The pressure is 
     * calculated for the box the integrator acts on and integrator's 
     * temperature is used for the ideal gas contribution.
     */
    public void setBox(Box newBox) {
        pc.setBox(newBox);
        box = newBox;
    }
    
    /**
     * Returns the integrator associated with this instance.  The pressure is 
     * calculated for the box the integrator acts on and integrator's 
     * temperature is used for the ideal gas contribution.
     */
    public Box getBox() {
        return box;
    }

    /**
     * Sets flag indicating whether calculated energy should include
     * long-range correction for potential truncation (true) or not (false).
     */
    public void setIncludeLrc(boolean b) {
    	iteratorDirective.includeLrc = b;
    }
    
    /**
     * Indicates whether calculated energy should include
     * long-range correction for potential truncation (true) or not (false).
     */
    public boolean isIncludeLrc() {
    	return iteratorDirective.includeLrc;
    }
    
    public void setTemperature(double newT) {
        pc.setTemperature(newT);
    }

	 /**
	  * Computes total pressure in box by summing virial over all pairs, and adding
	  * ideal-gas contribution.
	  */
    public IData getData() {
        if (box == null) {
            throw new IllegalStateException("You must call setBox before using this class");
        }
    	pc.zeroSum();
        potentialMaster.calculate(box, iteratorDirective, pc);
        data.x.E(pc.getPressureTensor());
        data.x.TE(1/box.getBoundary().volume());
        return data;
    }

    protected final DataTag tag;
    protected final DataTensor data;
    protected final DataInfo dataInfo;
    protected final PotentialMaster potentialMaster;
    protected Box box;
    protected IteratorDirective iteratorDirective;
    protected final PotentialCalculationPressureTensor pc;
    protected final double rD;
}
