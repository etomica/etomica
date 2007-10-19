package etomica.data.types;

import etomica.data.Data;
import etomica.data.DataTag;
import etomica.data.IDataInfo;
import etomica.data.IDataInfoFactory;
import etomica.units.Dimension;

public class DataPoint extends DataDoubleArray {

    private static final long serialVersionUID = 1L;

    public DataPoint(int[] arrayShape, double[] ptData) {
        super(arrayShape, ptData);
    }

    /**
     * Forms a DataFunction with the given array shape.
     */
    public DataPoint(int[] arrayShape) {
        super(arrayShape);
    }

    public static class DataInfoPoint extends DataInfoDoubleArray {

        private static final long serialVersionUID = 1L;

        public DataInfoPoint(String label, Dimension dimension, int[] arrayShape) {
            super(label, dimension, arrayShape);
        }
        
        public IDataInfoFactory getFactory() {
            return new DataInfoPointFactory(this);
        }
        
        public Data makeData() {
            return new DataFunction(arrayShape);
        }

    }
    
    public static class DataInfoPointFactory extends DataInfoDoubleArrayFactory {

        private static final long serialVersionUID = 1L;

        protected DataInfoPointFactory(DataInfoPoint template) {
            super(template);
        }
        
        public IDataInfo makeDataInfo() {
            DataInfoPoint dataInfo = new DataInfoPoint(label, dimension, getArrayShape());
            DataTag[] tagArray = new DataTag[tags.size()];
            dataInfo.addTags((DataTag[])tags.toArray(tagArray));
            return dataInfo;
        }

    }

}
