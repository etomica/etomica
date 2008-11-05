package etomica.normalmode;

import etomica.api.IData;
import etomica.data.DataPipe;
import etomica.data.DataProcessor;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.units.Null;

public class DataProcessorBoltzmannFactor extends DataProcessor {

	public DataProcessorBoltzmannFactor() {
		
		data = new DataDouble();
		dataInfo = new DataInfoDouble("Boltzmann Factor", Null.DIMENSION);
	}

	protected IData processData(IData inputData) {
		data.x = Math.exp(-inputData.getValue(0)/temperature);
		return data;
	}

	protected IEtomicaDataInfo processDataInfo(IEtomicaDataInfo inputDataInfo) {
		
		return dataInfo;
	}

	public DataPipe getDataCaster(IEtomicaDataInfo dataInfo) {
		return null;
	}
	
	public double getTemperature() {
		return temperature;
	}

	public void setTemperature(double temperature) {
		this.temperature = temperature;
	}

	protected final DataDouble data;  
	protected final DataInfoDouble dataInfo;
	protected double temperature;

}
