/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.action.BoxInflate;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.integrator.IntegratorBox;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Null;
import etomica.units.dimensions.Pressure;
import etomica.units.dimensions.Volume;

/**
 * Evaluates the pressure by examining the change in energy accompanying
 * small changes in volume.
 *
 * Pressure is approximately log(<x>)/deltaV
 * where x is one of the values returned by getData (<x> is that value averaged)
 * and deltaV is the change in volume associated with that x value.  deltaV = s * V
 * where s is the scaling factor (available via getScalingDataSource)
 */
public class MeterPressureByVolumeChange implements IDataSource {
    
    public MeterPressureByVolumeChange(Space space) {
        this(space, makeDefaultDimensions(space.D()));
    }
    
    private static final boolean[] makeDefaultDimensions(int D) {
        boolean[] dim = new boolean[D];
        for (int i=0; i<D; i++) {
            dim[i] = true;
        }
        return dim;
    }
    
    public MeterPressureByVolumeChange(Space space, boolean[] dimensions) {
        this(space, dimensions, new BoxInflate(space));
    }
    
    public MeterPressureByVolumeChange(Space space, BoxInflate pid){
        this(space, makeDefaultDimensions(space.D()), pid);
    }
    
    public MeterPressureByVolumeChange(Space space, boolean[] dimensions, BoxInflate pid){
        this.space = space;
        tag = new DataTag();
        setX(-0.001, 0.001, 10);
        inflateDimensions = new boolean[space.D()];
        setInflateDimensions(dimensions);
        inflater = pid;
        scale = space.makeVector();
    }

    /**
     * Sets the integrator associated with this instance.  The pressure is 
     * calculated for the box the integrator acts on and the integrator's 
     * temperature is used to calculate the pressure.
     */
    public void setIntegrator(IntegratorBox newIntegrator) {
        integrator = newIntegrator;
    }
    
    /**
     * Returns the integrator associated with this instance.  The pressure is 
     * calculated for the box the integrator acts on and the integrator's 
     * temperature is used to calculate the pressure.
     */
    public IntegratorBox getIntegrator() {
        return integrator;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }

    /**
     * For anisotropic volume change, indicates dimension in which volume is perturbed.
     */
    public final void setInflateDimensions(boolean[] directions) {
        if(directions.length != space.D()) {
            throw new IllegalArgumentException();
        }
        nDimension = 0;
        for(int i=0; i<directions.length; i++) {
            inflateDimensions[i] = directions[i];
            if(inflateDimensions[i]) nDimension++;
        }
        vDataSource.reset();
    }

    /**
     * Accessor method for setInflateDimension.
     */
    public boolean[] getInflateDimensions() {return inflateDimensions;}
    
    public void setX(double min, double max, int n) {
        if (n < 1) throw new IllegalArgumentException("n must be greater than 1");
        xDataSource = new DataSourceUniform("x", Volume.dimension(space.D()), n, min, max);
        vDataSource = new DataSourceExp(xDataSource);
        vDataSource.reset();
        dataInfo = new DataInfoDoubleArray("Pressure by Volume Change", Pressure.dimension(space.D()), new int[]{vDataSource.getDataInfo().getLength()});
        data = new DataDoubleArray(dataInfo.getLength());
        dataInfo.addTag(tag);
        dataArray = data.getData();
    }
    
    /**
     * Returns the data source for volume scalings.
     */
    public IDataSource getScalingDataSource() {
        return vDataSource;
    }
    
    public IData getData() {
        if (integrator == null) throw new IllegalStateException("must call setIntegrator before using meter");
        Box box = integrator.getBox();
        inflater.setBox(box);

        double uOld = integrator.getPotentialCompute().computeAll(false);
        final double[] x = ((DataDoubleArray)vDataSource.getData()).getData();
        double mult = 1.0/nDimension;
        for(int i=0; i<dataArray.length; i++) {
            scale.E(Math.pow(x[i],mult));
            for(int j=0; j<space.D(); j++) {
                if(!inflateDimensions[j]) scale.setX(j,1.0);
            }

            inflater.setVectorScale(scale);
            inflater.actionPerformed();
            double uNew = integrator.getPotentialCompute().computeAll(false);
            dataArray[i] = Math.exp(-(uNew-uOld)/integrator.getTemperature()
                              + box.getMoleculeList().size()*(x[i]-1));
            inflater.undo();
        }

        return data;
    }

    private DataDoubleArray data;
    private IDataInfo dataInfo;
    private final DataTag tag;
    private double[] dataArray;
    private final BoxInflate inflater;
    private final Vector scale;
    private final boolean[] inflateDimensions;
    private int nDimension;
    private final Space space;
    private DataSourceUniform xDataSource;
    protected DataSourceExp vDataSource;
    private IntegratorBox integrator;
    
    /**
     * Transforms the scaling from linear (-s to +s) to exponential (exp(-s) to exp(+s))
     */
    protected static class DataSourceExp implements IDataSource {
        public DataSourceExp(DataSourceUniform wrappedDataSource) {
            this.wrappedDataSource = wrappedDataSource;
            tag = new DataTag();
            int n = wrappedDataSource.getDataInfo().getLength();
            if (n % 2 == 1) n--;
            data = new DataDoubleArray(n);
            dataInfo = new DataInfoDoubleArray("something", Null.DIMENSION, new int[]{n});
        }
        
        public DataTag getTag() {
            return tag;
        }
        
        public IDataInfo getDataInfo() {
            return dataInfo;
        }
        
        public void reset() {
            double[] x = ((DataDoubleArray)wrappedDataSource.getData()).getData();
            double[] scaledX = data.getData();
            int xIndex = 0;
            double dx = (wrappedDataSource.getXMax()-wrappedDataSource.getXMin())/wrappedDataSource.getNValues();
            for (int i=0; i<scaledX.length; i++, xIndex++) {
                // skip x=0 in the center
                if (Math.abs(x[xIndex]) < 0.25*dx) xIndex++;
                scaledX[i] = Math.exp(x[xIndex]);
            }
        }
        
        public IData getData() {
            return data;
        }
        
        protected final DataSourceUniform wrappedDataSource;
        protected final DataTag tag;
        protected DataInfoDoubleArray dataInfo;
        protected DataDoubleArray data;
    }
}
