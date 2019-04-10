package etomica.cavity;

import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataGroup;
import etomica.math.numerical.PolynomialFit;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Null;

class DataProcessorFit extends DataProcessorForked {
    protected String label;
    protected final int order;
    protected final boolean log;
    protected final double[] x, y, w;
    protected final DataFunction data;
    protected DataDoubleArray xData;
    protected double xMin, xMax;
    protected double[] eDataSave;

    public DataProcessorFit(String label, int nPoints, int order, boolean log) {
        this(label, nPoints, order, log, 0, 1);
    }

    public DataProcessorFit(String label, int nPoints, int order, boolean log, double xMin, double xMax) {
        this.label = label;
        this.order = order;
        this.log = log;
        this.xMin = xMin;
        this.xMax = xMax;
        x = new double[nPoints];
        y = new double[nPoints];
        w = new double[nPoints];
        DataSourceUniform xDataSource = new DataSourceUniform();
        xDataSource.setTypeMin(DataSourceUniform.LimitType.INCLUSIVE);
        xDataSource.setTypeMax(DataSourceUniform.LimitType.INCLUSIVE);
        xDataSource.setNValues(101);
        xDataSource.setXMin(xMin);
        xDataSource.setXMax(xMax);
        data = new DataFunction(new int[]{xDataSource.getNValues()});
        double[] xOut = ((DataDoubleArray) xDataSource.getData()).getData();
        DataSourceIndependentSimple rData = new DataSourceIndependentSimple(xOut, new DataDoubleArray.DataInfoDoubleArray("r", Length.DIMENSION, new int[]{xOut.length}));
        dataInfo = new DataFunction.DataInfoFunction(label, Null.DIMENSION, rData);
        dataInfo.addTag(tag);
    }

    @Override
    protected IData processData(IData inputData) {
        IData yData = ((DataGroup) inputData).getData(0);
        IData eData = ((DataGroup) inputData).getData(1);

        int i = 0;
        int nGood = 0;
        for (int j = 0; i < x.length && j < yData.getLength(); j++) {
            if (yData.getValue(j) == 0) continue;
            x[i] = xData.getValue(j);
            if (x[i] < xMin || x[i] > xMax) continue;
//            System.out.println(x[i]+" "+y[i]+" "+w[i]);
            if (log) {
                y[i] = Math.log(yData.getValue(j));
                double ratio = eData.getValue(j) / yData.getValue(j);
                w[i] = 1 / (ratio * ratio);
                if (ratio > 0.2) w[i] = 0;
            } else {
                y[i] = yData.getValue(j);
                w[i] = 1.0 / (eData.getValue(j) * eData.getValue(j));
            }
            if (Double.isNaN(w[i]) || w[i] == 0) w[i] = 0;
            else nGood++;
            i++;
        }
        double[] yOut = data.getData();
//            System.out.println("Found "+nGood+" good points");
        if (nGood < 4) {
            eDataSave = new double[yOut.length];
            for (int j = 0; j < yOut.length; j++) yOut[j] = eDataSave[j] = Double.NaN;
            return data;
        }
        for (; i < x.length; i++) {
            x[i] = y[i] = w[i] = 0;
        }
        PolynomialFit.FitResult fr = null;
        for (int o = 1; o <= order && o < nGood * 2; o++) {
            fr = PolynomialFit.doFit(o, x, y, w, true);
            double[] poly = fr.coeff;
            double chi = PolynomialFit.getChi(x, y, w, poly);
//            System.out.println(o+" chi " + chi);
            if (chi < 1) break;
        }
        DataDoubleArray rData = ((DataFunction.DataInfoFunction) dataInfo).getXDataSource().getIndependentData(0);
        double[][] yFit = PolynomialFit.getFit(rData.getData(), fr);
        eDataSave = new double[yFit[1].length];
        for (int j = 0; j < data.getLength(); j++) {
            yOut[j] = yFit[0][j];
            eDataSave[j] = yFit[1][j];
            if (log) {
                eDataSave[j] /= yOut[j];
                yOut[j] = Math.exp(yOut[j]);
            }
        }
        return data;
    }

    public double[] getLastErr() {
        return eDataSave;
    }

    @Override
    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        xData = ((DataFunction.DataInfoFunction) (((DataGroup.DataInfoGroup) inputDataInfo).getSubDataInfo(0))).getXDataSource().getIndependentData(0);
        return dataInfo;
    }
}
