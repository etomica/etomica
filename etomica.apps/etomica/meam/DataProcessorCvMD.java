/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.meam;

import etomica.data.DataPipe;
import etomica.data.DataProcessor;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.units.CompoundDimension;
import etomica.units.Dimension;
import etomica.units.Energy;
import etomica.units.Quantity;
import etomica.units.Temperature;


public class DataProcessorCvMD extends DataProcessor {

    private static final long serialVersionUID = 1L;

    IntegratorVelocityVerlet integrator;

	DataDouble data;

	DataInfoDouble dataInfo;

	public void setIntegrator(IntegratorVelocityVerlet integrator) {
		this.integrator = integrator;
	}

	public IData processData(IData inputData){
		data.x = ((DataDouble)inputData).x;
		double systemTemp = integrator.getTemperature();
		data.x /= systemTemp;
		data.x *= data.x/integrator.getBox().getMoleculeList().getMoleculeCount();
		return data;
	}

	public DataPipe getDataCaster(IEtomicaDataInfo inputDataInfo){
		return null;
	}

	public IEtomicaDataInfo processDataInfo(IEtomicaDataInfo inputDataInfo){
		Dimension cvDimension = new CompoundDimension(new Dimension[]{Energy.DIMENSION, Temperature.DIMENSION, Quantity.DIMENSION}, new double[]{1,-1,-1});
		data = new DataDouble();
		dataInfo = new DataInfoDouble("Heat Capacity", cvDimension);
		return dataInfo;
	}

}