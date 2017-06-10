/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.densityofstate;

import etomica.data.DataPipe;
import etomica.data.DataProcessor;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.units.dimensions.Null;

public class DataProcessorPhi extends DataProcessor {
	
	public DataProcessorPhi(){
		tag = new DataTag();
	}

	protected IData processData(IData inputData) {
		// TODO Auto-generated method stub
		double U = ((DataDouble)inputData).x;
		double [] phi = data.getData();
		
		double UT = U/temperature;
		double UT2 = UT*UT;
		double UT4 = UT2*UT2;
		double UT6 = UT4*UT2;
		double UT8 = UT4*UT4;
		
		phi[0] = 1/Math.sqrt(temperature);
		phi[1] =(1-  UT)/Math.sqrt(temperature);
		phi[2] =(1-2*UT+ (1/2)*UT2)/Math.sqrt(temperature);
		phi[3] =(1-3*UT+ (3/2)*UT2- (1/6)*UT*UT2)/Math.sqrt(temperature);
		phi[4] =(1-4*UT+   3  *UT2- (2/3)*UT*UT2+ (1/24)*UT4)/Math.sqrt(temperature);
		phi[5] =(1-5*UT+   5  *UT2- (5/3)*UT*UT2+ (5/24)*UT4-(1/120)*UT*UT4)/Math.sqrt(temperature);
		phi[6] =(1-6*UT+(15/2)*UT2-(10/3)*UT*UT2+ (5/8) *UT4- (1/20)*UT*UT4+(1/720)*UT6)/Math.sqrt(temperature);
		phi[7] =(1-7*UT+(21/2)*UT2-(35/6)*UT*UT2+(35/24)*UT4- (7/40)*UT*UT4+(7/720)*UT6-(1/5040)*UT*UT6)/Math.sqrt(temperature);
		phi[8] =(1-8*UT+  14  *UT2-(28/3)*UT*UT2+(35/12)*UT4- (7/15)*UT*UT4+(7/180)*UT6- (1/630)*UT*UT6+(1/40320)*UT8)/Math.sqrt(temperature);
		phi[9] =(1-9*UT+  18  *UT2-  14  *UT*UT2+ (21/4)*UT4-(21/20)*UT*UT4+( 7/60)*UT6- (1/140)*UT*UT6+ (1/4480)*UT8-(1/362880)*UT*UT8)/Math.sqrt(temperature);
		phi[10] = phi[0]*phi[0];
		phi[11] = phi[1]*phi[1];
		phi[12] = phi[2]*phi[2];
		phi[13] = phi[3]*phi[3];
		phi[14] = phi[4]*phi[4];
		phi[15] = phi[5]*phi[5];
		phi[16] = phi[6]*phi[6];
		phi[17] = phi[7]*phi[7];
		phi[18] = phi[8]*phi[8];
		phi[19] = phi[9]*phi[9];
		
		return data;
		
	}

	public void setTemperature(double newTemperature){
		temperature = newTemperature;
	}
	protected IEtomicaDataInfo processDataInfo(IEtomicaDataInfo inputDataInfo) {
		// TODO Auto-generated method stub
		dataInfo = new DataInfoDoubleArray("phi", Null.DIMENSION, new int []{20});
		data = new DataDoubleArray(20);
		return dataInfo;
	}

	public DataPipe getDataCaster(IEtomicaDataInfo dataInfo) { //hook up to meter to the Dataprocessor
		if (dataInfo instanceof DataInfoDouble)
			return null;
		throw new IllegalArgumentException("i only want double");
	}
	
	private DataTag tag;
	private DataInfoDoubleArray dataInfo;
	private DataDoubleArray data;
	private double temperature;
	
}
