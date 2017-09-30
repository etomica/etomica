/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.types;

import com.fasterxml.jackson.annotation.JsonIgnore;
import com.fasterxml.jackson.annotation.JsonProperty;
import etomica.data.DataSourceIndependent;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataInfoFactory;
import etomica.units.dimensions.Dimension;


/**
 * Collects two or more DataDoubleArray instances, and organizes them into
 * "dependent" and "independent" subgroups. Normally these data represent a
 * functional dependence, in which the dependent data were calculated at each
 * point in the domain of the independent data. However, nothing about this
 * class enforces the dependence, it merely classifies the data into the two
 * groups.
 * <p>
 * Multidimensional functions hold the independent data as an array of
 * one-dimensional DataDoubleArrays, and the dependent data is represented by a
 * single DataDoubleArray that is formed with a corresponding number of
 * dimensions.
 * <p>
 * All arithmetic operations apply only to the dependent data. Independent data
 * values are unaffected by them.
 * <p>
 * Note that all instances created by the same factory will hold the same instances
 * of the independent data.  Thus if the independent data is changed in any instance,
 * it will be reflected in all instances that were constructed by the factory.
 * 
 * @author David Kofke and Andrew Schultz
 *  
 */

public class DataFunction extends DataDoubleArray {

    /**
     * Forms a DataFunction using the given independent and dependent data
     * instances.
     *
     * @param yData
     *            the data defined on the space of independent data; a
     *            double[] with getDimension equal to
     *            independentData.length
     * 
     * @throws IllegalArgumentException
     *             if length of independentData array does not equal
     *             dependentData.getArrayDimension(), or if any of the
     *             independent data have array dimension not equal to 1
     */
    public DataFunction(int[] arrayShape, double[] yData) {
        super(arrayShape, yData);
    }
    
    /**
     * Forms a DataFunction with the given array shape.
     */
    public DataFunction(int[] arrayShape) {
        super(arrayShape);
    }

    private static final long serialVersionUID = 1L;

    public static class DataInfoFunction extends DataInfoDoubleArray {
        public DataInfoFunction(String label, Dimension dimension, DataSourceIndependent xDataSource) {
            super(label, dimension, getArrayShape(xDataSource));
            this.xDataSource = xDataSource;
        }
        
        private static int[] getArrayShape(DataSourceIndependent xDataSource) {
            int[] arrayShape = new int[xDataSource.getIndependentArrayDimension()];
            for (int i=0; i<arrayShape.length; i++) {
                arrayShape[i] = xDataSource.getIndependentDataInfo(i).getArrayShape()[0];
            }
            return arrayShape;
        }

        @JsonIgnore
        public DataSourceIndependent getXDataSource() {
            return xDataSource;
        }

        @JsonIgnore
        public IEtomicaDataInfoFactory getFactory() {
            return new DataInfoFunctionFactory(this);
        }
        
        public IData makeData() {
            return new DataFunction(arrayShape);
        }

        protected final DataSourceIndependent xDataSource;

        @JsonProperty
        private double[] getIndependentData() {
            return getXDataSource().getIndependentData(0).getData();
        }

        @JsonProperty
        private IEtomicaDataInfo getIndependentDataInfo() {
            return getXDataSource().getIndependentDataInfo(0);
        }
    }
    
    public static class DataInfoFunctionFactory extends DataInfoDoubleArrayFactory {
        protected DataInfoFunctionFactory(DataInfoFunction template) {
            super(template);
            xDataSource = template.xDataSource;
        }
        
        public IEtomicaDataInfo makeDataInfo() {
            DataInfoFunction dataInfo = new DataInfoFunction(label, dimension, xDataSource);
            DataTag[] tagArray = new DataTag[tags.size()];
            dataInfo.addTags((DataTag[])tags.toArray(tagArray));
            return dataInfo;
        }
        
        /**
         * Sets array of independent DataInfo.  The array is copied so further
         * changes made to the given array will not be affect this factory.
         */
        public void setXDataSource(DataSourceIndependent newXDataSource) {
            xDataSource = newXDataSource;
        }
        
        /**
         * Returns a copy of array of independent DataInfo.
         */
        public DataSourceIndependent getXDataSource() {
            return xDataSource;
        }
        
        private static final long serialVersionUID = 1L;
        protected DataSourceIndependent xDataSource;
    }
}
