package etomica.models.oneDHardRods;

import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataSource;
import etomica.data.DataTag;
import etomica.data.IDataInfo;

public class MeterNMConversion implements DataSource {

	
	Data data;
	DataInfo dataInfo;
	DataTag tag;
	int nTests;		//the number of tests made
	
	
	public MeterNMConversion(){
		
		
		
	}
	
	
	
	
	private void insertTest(){
		
	}
	
	private void removeTest(){
		
	}
	
	
	public Data getData(){
		

		//bits that choose insertion or removal of a mode
		//probably some sort of degrees of freedom check here
		
		//if else statement to call insertTest or removeTest
		// inside the if statement for loop to repeatedly test?
		
		
		
		return data;
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
