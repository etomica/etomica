package etomica.data.types;

import java.io.Serializable;

import etomica.Data;
import etomica.DataInfo;
import etomica.data.DataFactory;
import etomica.data.DataProcessor;
import etomica.space.Tensor;

/*
 * History
 * Created on Jul 28, 2005 by kofke
 */
public class CastToTable extends DataProcessor implements Serializable {

        public DataInfo processDataInfo(DataInfo inputDataInfo) {
            Class inputClass = inputDataInfo.getClass();
            DataFactory factory = inputDataInfo.getDataFactory();
//            DataTable.Column[] myColumns = null;
            if (inputClass == DataDoubleArray.class) {
                inputType = 0;
                int[] arrayShape = ((DataDoubleArray.Factory)factory).getArrayShape();
                if(arrayShape.length > 2) {
                    throw new IllegalArgumentException("Cannot cast to table a data set with dimension greater than 2: "+arrayShape.length);
                }
                if(arrayShape.length == 1) {
                    outputData = new DataTable(inputDataInfo.getLabel(),new DataInfo[]{inputDataInfo},arrayShape[0]);
                } else {
                    outputData = new DataTable(inputDataInfo.getLabel(),inputDataInfo.getDimension(),arrayShape[0],arrayShape[1]);
                }
            } else if (inputClass == DataDouble.class) {
                inputType = 1;
                outputData = new DataTable(inputDataInfo.getLabel(),new DataInfo[]{inputDataInfo},1);
            } else if (inputClass == DataInteger.class) {
                inputType = 2;
                outputData = new DataTable(inputDataInfo.getLabel(),new DataInfo[]{inputDataInfo},1);
            } else if (inputClass == DataVector.class) {
                inputType = 3;
                int D = ((DataVector.Factory)factory).getSpace().D();
                outputData = new DataTable(inputDataInfo.getLabel(),new DataInfo[]{inputDataInfo},D);
            } else if (inputClass == DataTensor.class) {
                inputType = 4;
                int D = ((DataTensor.Factory)factory).getSpace().D();
                outputData = new DataTable(inputDataInfo.getLabel(),inputDataInfo.getDimension(),D,D);
            } else if(inputClass == DataFunction.class) {
                int[] sizes = ((DataFunction.Factory)factory).getIndependentDataSizes();
                DataDoubleArray[] data = ((DataFunction.Factory)factory).getIndependentData();
                if (sizes.length == 1) {
                    outputData = new DataTable(inputDataInfo.getLabel(),new DataInfo[]{data[0].getDataInfo(),inputDataInfo},sizes[0]);
                }
                else {//if (sizes.length == 2) {
                    throw new IllegalArgumentException("DataFunction must be 1-dimensional");
                }
                inputType = 5;
            } else {
                throw new IllegalArgumentException("Cannot cast to DataTable from "
                        + inputClass);
            }
            return outputData.getDataInfo();

        }//end of processDataInfo
        
        protected Data processData(Data data) {
            switch(inputType) {
            case 0:
                int nColumns = outputData.myColumns.length;
                for (int i=0; i<nColumns; i++) {
                    ((DataDoubleArray)data).assignColumnTo(i,outputData.myColumns[i].data);
                }
                break;
            case 1:
                outputData.myColumns[0].data[0] = ((DataDouble)data).x;
                break;
            case 2:
                outputData.myColumns[0].data[0] = ((DataInteger)data).x;
                break;
            case 3:
                ((DataVector)data).x.assignTo(outputData.myColumns[0].data);
                break;
            case 4:
                Tensor x = ((DataTensor)data).x;
                for (int i=0; i<x.length(); i++) {
                    for (int j=0; j<x.length(); j++) {
                        outputData.myColumns[i].data[j] = x.component(i,j);
                    }
                }
                break;
            case 5:
                ((DataFunction)data).getXData(0).assignColumnTo(0,outputData.myColumns[0].data);
                ((DataFunction)data).getYData().assignColumnTo(0,outputData.myColumns[1].data);
                break;
            }
            return outputData;
        }
        
        public DataProcessor getDataCaster(DataInfo dataInfo) {
            return null;
        }
                
        DataTable outputData;
        DataDoubleArray arrayData;
        int inputType;
    }