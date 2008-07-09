package etomica.normalmode;

import etomica.data.Data;
import etomica.data.DataPipe;
import etomica.data.DataProcessor;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.units.Null;

public class DataProcessorBoltzmannFactor extends DataProcessor {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	public DataProcessorBoltzmannFactor() {
		
		data = new DataDouble();
		dataInfo = new DataInfoDouble("Boltzmann Factor", Null.DIMENSION);
	}

	protected Data processData(Data inputData) {
		data.x = Math.exp(-inputData.getValue(0)/temperature);
		return data;
	}

	protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
		
		return dataInfo;
	}

	public DataPipe getDataCaster(IDataInfo dataInfo) {
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
