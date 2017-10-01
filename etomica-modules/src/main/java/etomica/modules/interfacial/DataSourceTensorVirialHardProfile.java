/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.interfacial;

import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.data.types.DataGroup;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.integrator.IntegratorHard;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Length;

import java.util.Arrays;

/**
 * A DataSource that collects the average virial tensor as a function of
 * position for one dimension.  The profile data is exposed via an inner class
 * that fronts that data (since this class extends DataSourceTensorVirialHard).
 * 
 * @author Andrew Schultz
 */
public class DataSourceTensorVirialHardProfile extends DataSourceTensorVirialHard implements DataSourceIndependent {
    
    public DataSourceTensorVirialHardProfile(Space space) {
        super(space);
        profileData = new DataGroup(new IData[0]);
        profileDataInfo = new DataInfoGroup("Virial profiles", Energy.DIMENSION, new IEtomicaDataInfo[0]);
        profileDataInfo.addTag(tag);
        binSize = 0.1;
        virialProfile = new double[space.D()][0];
        xTag = new DataTag();
    }
    
    /**
     * Current value of the meter, obtained by dividing sum of collision virial contributions by time elapsed since last call.
     * If elapsed-time interval is zero, returns the value reported at the last call to the method.
     */
    public DataGroup getVirialProfile() {
        double currentTime = integratorHard.getCurrentTime();
        double elapsedTime = currentTime - lastProfileTime;
        lastProfileTime = currentTime;
        
        Vector boxDim = integratorHard.getBox().getBoundary().getBoxSize();
        if (L != boxDim.getX(0)) {
            // the data we collected is bogus.  reset and return NaN.
            setupProfileData();
            profileData.E(Double.NaN);
            return profileData;
        }
        
        if(elapsedTime == 0.0) {
            // oops
            profileData.E(Double.NaN);
            return profileData;
        }

        for (int i=0; i<virialProfile.length; i++) {
            System.arraycopy(virialProfileWork[i], 0, virialProfile[i], 0, virialProfile[i].length);
            Arrays.fill(virialProfileWork[i], 0);
        }
        profileData.TE(-1.0/elapsedTime);
        return profileData;
    }
    
    public DataInfoGroup getProfileDataInfo() {
        return profileDataInfo;
    }
    
    /**
     * Sums contribution to virial for each collision.
     */
    public void collisionAction(IntegratorHard.Agent agent) {
        // superclass collects total virial
        super.collisionAction(agent);

        Tensor virialTensor = agent.collisionPotential.lastCollisionVirialTensor();

        double x0 = agent.atom.getPosition().getX(0);
        // wrap around PBC
        x0 -= Math.round(x0*Li) * L;
        int iBin0 = (int) ((x0 + halfL) / binSize);
        double x1 = agent.collisionPartner.getPosition().getX(0);
        // wrap around PBC
        x1 -= Math.round(x1*Li) * L;
        int iBin1 = (int) ((x1 + halfL) / binSize);
        if (iBin0 == iBin1) {
            // both in 1 bin
            addTensor(iBin0, virialTensor);
        }
        else if (Math.abs(iBin1-iBin0) > nBins/2) {
            // wrapped around the box
            if (iBin0 > nBins/2) {
                //swap so iBin1 is on the right side, iBin0 on the left
                int foo = iBin0;
                iBin0 = iBin1;
                iBin1 = foo;
            }
            // bins on the right: (nBins-1)-iBin1+1
            // bins on the left:  iBin0+1
            int nBetween = (nBins-1)-iBin1+1 + iBin0+1;
            virialTensor.TE(1.0/nBetween);
            for (int i=iBin1; i<nBins; i++) {
                addTensor(i, virialTensor);
            }
            for (int i=0; i<iBin0+1; i++) {
                addTensor(i, virialTensor);
            }
        }
        else {
            if (iBin0 > iBin1) {
                //swap so iBin1 is on the right side, iBin0 on the left
                int foo = iBin0;
                iBin0 = iBin1;
                iBin1 = foo;
            }
            int nBetween = iBin1 - iBin0 + 1;
            virialTensor.TE(1.0/nBetween);
            for (int i=iBin0; i<iBin1+1; i++) {
                addTensor(i, virialTensor);
            }
        }
    }

    /**
     * Adds the given contribution to the virial sum for bin iBin.
     */
    protected void addTensor(int iBin, Tensor virialTensor) {
        for (int i=0; i<virialTensor.D(); i++) {
            virialProfileWork[i][iBin] += virialTensor.component(i,i);
        }
    }
    
    public void setBinSize(double newBinSize) {
        binSize = newBinSize;
        if (integratorHard != null) {
            setupProfileData();
        }
    }
    
    public double getBinSize() {
        return binSize;
    }
    
    public void setIntegrator(IntegratorHard newIntegrator) {
        super.setIntegrator(newIntegrator);
        setupProfileData();
    }

    protected void setupProfileData() {
        Vector boxDim = integratorHard.getBox().getBoundary().getBoxSize();
        L = boxDim.getX(0);
        Li = 1.0/L;
        halfL = 0.5*L;
        nBins = (int)Math.round(L/binSize);
        binSize = boxDim.getX(0) / nBins;
        
        xData = new DataDoubleArray(nBins);
        double[] x = xData.getData();
        for (int i=0; i<nBins; i++) {
            x[i] = -0.5 * boxDim.getX(0) + (i+0.5) * binSize;
        }
        xDataInfo = new DataInfoDoubleArray("x", Length.DIMENSION, new int[] {nBins});
        
        DataFunction[] virialData = new DataFunction[boxDim.getD()];
        DataInfoFunction[] virialDataInfo = new DataInfoFunction[boxDim.getD()];
        for (int i=0; i<virialData.length; i++) {
            virialData[i] = new DataFunction(new int[]{nBins});
            virialProfile[i] = virialData[i].getData();
            virialDataInfo[i] = new DataInfoFunction("virial", Energy.DIMENSION, this);
        }
        profileData = new DataGroup(virialData);
        profileDataInfo = new DataInfoGroup("virial", Energy.DIMENSION, virialDataInfo);
        
        virialProfileWork = new double[boxDim.getD()][nBins];
    }

    public int getIndependentArrayDimension() {
        return 1;
    }

    public DataDoubleArray getIndependentData(int i) {
        return xData;
    }

    public DataInfoDoubleArray getIndependentDataInfo(int i) {
        return xDataInfo;
    }

    public DataTag getIndependentTag() {
        return xTag;
    }

    private static final long serialVersionUID = 1L;
    protected double lastProfileTime;
    protected DataGroup profileData;
    protected double[][] virialProfile, virialProfileWork;
    protected DataInfoGroup profileDataInfo;
    protected DataDoubleArray xData;
    protected DataInfoDoubleArray xDataInfo;
    protected final DataTag xTag;
    protected double binSize;
    protected double L, Li, halfL;
    protected int nBins;

    public static class DataSourceVirialProfile implements IDataSource {

        public DataSourceVirialProfile(DataSourceTensorVirialHardProfile meter) {
            this.meter = meter;
            tag = new DataTag();
        }
        
        public IData getData() {
            return meter.getVirialProfile();
        }

        public IEtomicaDataInfo getDataInfo() {
            return meter.getProfileDataInfo();
        }

        public DataTag getTag() {
            return tag;
        }
        
        protected final DataSourceTensorVirialHardProfile meter;
        protected final DataTag tag;
    }
}
