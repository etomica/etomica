/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.interfacial;

import etomica.box.Box;
import etomica.space.Vector;
import etomica.data.DataPipe;
import etomica.data.DataProcessor;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataInfoFactory;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.data.types.DataGroup;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.space.Space;
import etomica.units.dimensions.Area;
import etomica.units.dimensions.DimensionRatio;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Length;

public class DataProcessorInterfacialTensionProfile extends DataProcessor {

    public DataProcessorInterfacialTensionProfile(Space space) {
        virialData = new double[space.D()][0];
    }
    
    public void setBox(Box newBox) {
        box = newBox;
    }

    public Box getBox() {
        return box;
    }

    public void setProfileDim(int newProfileDim) {
        profileDim = newProfileDim;
    }

    public int getProfileDim() {
        return profileDim;
    }

    protected IData processData(IData inputData) {
        DataGroup dataGroup = (DataGroup)inputData;
        int D = virialData.length;
        for (int i=0; i<D; i++) {
            virialData[i] = ((DataFunction)dataGroup.getData(i)).getData();
        }
        int nBins = data.getArrayShape(0);
        double[] tension = data.getData();
        for (int i=0; i<nBins; i++) {
            tension[i] = (D-1)*virialData[profileDim][i];
        }
        for (int j=0; j<D; j++) {
            if (j == profileDim) continue;
            for (int i=0; i<nBins; i++) {
                tension[i] -= virialData[j][i];
            }
        }

        double area = 1;
        Vector dim = box.getBoundary().getBoxSize();
        for (int i=0; i<dim.getD(); i++) {
            if (i == profileDim) continue;
            area *= dim.getX(i);
        }
        double binSize = dim.getX(0) / virialData[profileDim].length;
        double fac = 0.5/area/binSize/(D-1);
        for (int i=0; i<nBins; i++) {
            tension[i] *= fac;
        }
        return data;
    }

    protected IEtomicaDataInfo processDataInfo(IEtomicaDataInfo inputDataInfo) {
        DataInfoFunction dataInfo0 = (DataInfoFunction)((DataInfoGroup)inputDataInfo).getSubDataInfo(0);
        data = (DataFunction)dataInfo0.makeData();
        IEtomicaDataInfoFactory dataInfoFactory = dataInfo0.getFactory();
        dataInfoFactory.setDimension(new DimensionRatio(Energy.DIMENSION, ((DataInfoGroup)inputDataInfo).getNDataInfo() == 2 ? Length.DIMENSION : Area.DIMENSION));
        dataInfoFactory.setLabel("Interfacial tension profile");
        return dataInfoFactory.makeDataInfo();
    }

    public DataPipe getDataCaster(IEtomicaDataInfo inputDataInfo) {
        return null;
    }

    protected DataFunction data;
    protected final double[][] virialData;
    protected Box box;
    protected int profileDim;
}
