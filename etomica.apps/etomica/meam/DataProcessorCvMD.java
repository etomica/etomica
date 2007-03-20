/*
 * Created on Aug 8, 2006
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.meam;

import etomica.data.Data;
import etomica.data.DataPipe;
import etomica.data.DataProcessor;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.units.CompoundDimension;
import etomica.units.Dimension;
import etomica.units.Energy;
import etomica.units.Quantity;
import etomica.units.Temperature;


final class DataProcessorCvMD extends DataProcessor {
	IntegratorVelocityVerlet integrator;

	DataDouble data;

	DataInfoDouble dataInfo;

	public void setIntegrator(IntegratorVelocityVerlet integrator) {
		this.integrator = integrator;
	}

	public Data processData(Data inputData){
		data.x = ((DataDouble)inputData).x;
		double systemTemp = integrator.getTemperature();
		data.x /= systemTemp;
		data.x *= data.x/integrator.getPhase().moleculeCount();
		return data;
	}

	public DataPipe getDataCaster(IDataInfo inputDataInfo){
		return null;
	}

	public IDataInfo processDataInfo(IDataInfo inputDataInfo){
		Dimension cvDimension = new CompoundDimension(new Dimension[]{Energy.DIMENSION, Temperature.DIMENSION, Quantity.DIMENSION}, new double[]{1,-1,-1});
		data = new DataDouble();
		dataInfo = new DataInfoDouble("Heat Capacity", cvDimension);
		return dataInfo;
	}
}