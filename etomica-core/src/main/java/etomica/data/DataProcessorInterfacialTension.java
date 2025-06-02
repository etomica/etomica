/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.box.Box;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Area;
import etomica.units.dimensions.DimensionRatio;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Length;

/**
 * Data Processor that takes the virial components as input data and returns
 * the interfacial tension.  This class must be told which dimension the
 * surfaces exist in.
 *
 * @author Andrew Schultz
 */
public class DataProcessorInterfacialTension extends DataProcessor {

    protected final Space space;
    protected final DataDouble data;
    protected Box box;
    protected int surfaceDim;

    public DataProcessorInterfacialTension(Space space) {
        this.space = space;
        data = new DataDouble();
    }

    public Box getBox() {
        return box;
    }

    public void setBox(Box newBox) {
        box = newBox;
    }

    /**
     * Returns the dimension of the surfaces.  0 (x) means that the surfaces
     * exist in the yz plane.
     */
    public int getSurfaceDimension() {
        return surfaceDim;
    }

    /**
     * Sets the dimension of the surfaces.  0 (x) means that the surfaces exist
     * in the yz plane.
     */
    public void setSurfaceDimension(int newDim) {
        surfaceDim = newDim;
    }

    protected IData processData(IData inputData) {
        double area = 1;
        Vector dim = box.getBoundary().getBoxSize();
        int D = dim.getD();
        data.x = inputData.getValue(surfaceDim);
        for (int i=0; i<D; i++) {
            if (i == surfaceDim) continue;
            area *= dim.getX(i);
            data.x -= inputData.getValue(i)/(D-1);
        }
        data.x /= 2*area;
        return data;
    }

    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        return new DataInfoDouble("Interfacial tension", new DimensionRatio(Energy.DIMENSION, space.D() == 2 ? Length.DIMENSION : Area.DIMENSION));
    }
}
