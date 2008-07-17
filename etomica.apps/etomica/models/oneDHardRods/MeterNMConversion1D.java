package etomica.models.oneDHardRods;

import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataSourceScalar;
import etomica.data.DataTag;
import etomica.data.IDataInfo;
import etomica.space.ISpace;
import etomica.units.Null;

public class MeterNMConversion1D extends DataSourceScalar {

	
	Data data;
	DataInfo dataInfo;
	DataTag tag;
	int nTests;		//the number of tests made
	
	
	public MeterNMConversion1D(ISpace space){
		super("exp(-\u03BC/kT)", Null.DIMENSION);//"\u03BC" is Unicode for greek "mu"
		nTests = 100;
		
		
	}
	
	
	
	
	private void insertTest(){
		
	}
	
	private void removeTest(){
		
	}
	
	
	public double getDataAsScalar(){
		double uTest = 0.0;

		//bits that choose insertion or removal of a mode
		//probably some sort of degrees of freedom check here
		
		//if else statement to call insertTest or removeTest
		// inside the if statement for loop to repeatedly test?
		
		
		return uTest;
	}
	
	public void setModes(int[] m){
		
	}
	
	 public IDataInfo getDataInfo(){
		 return (IDataInfo)dataInfo;
	 }
	 
	 public DataTag getTag(){
		 return tag;
	 }
	 public void setNTests(int n){
		 nTests = n;
	 }
	 public int getNTests(){
		 return nTests;
	 }
}
