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
	
	
	
	
	
	
	
	public Data getData(){
		
		
		
		
		return data;
	}
	
	 public IDataInfo getDataInfo(){
		 return (IDataInfo)dataInfo;
	 }
	 
	 public DataTag getTag(){
		 return tag;
	 }
}
