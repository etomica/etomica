package etomica.data;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction.DataInfoFunction;

/**
 * DataProcessor that interpolates incoming data using an Akima spline
 * interpolation scheme.
 * 
 *   "A New Method of Interpolation and Smooth Curve Fitting Based on Local
 *   Procedures", H. Akima, JACM 17 (1970), 589
 * 
 * @author Andrew Schultz
 */
public class AkimaSpline extends DataProcessor implements DataSourceIndependent {

    public AkimaSpline(int nSubPoints) {
        x = new double[0];
        y = new double[0];
        m = new double[0];
        t = new double[0];
        interpolatedX = new double[0];
        interpolatedY = new double[0];
        this.nSubPoints = nSubPoints;
        xTag = new DataTag();
    }
    
    protected void doInterpolation(double[] originalX, double[] originalY) {
        int N = originalX.length+4;
        if (x.length != N) {
            x = new double[N];
            y = new double[N];
            m = new double[N];
            t = new double[N];
        }
        System.arraycopy(originalX, 0, x, 2, N-4);
        System.arraycopy(originalY, 0, y, 2, N-4);

        for (int i=2; i<N-5; i++) {
            m[i] = (y[i+1]-y[i])/(x[i+1]-x[i]);
        }

        m[1] = 2*m[2] - m[3];
        m[0] = 2*m[1] - m[2];

        m[N-3] = 2*m[N-4] - m[N-5];
        m[N-2] = 2*m[N-3] - m[N-4];

        for (int i=2; i<N-2; i++) {
            double am2m1 = Math.abs(m[i-1]-m[i-2]);
            double am4m3 = Math.abs(m[i+1]-m[i]);
            if (am2m1 != 0 || am4m3 != 0 || am2m1 == am4m3) {
                t[i] = (am4m3*m[i-1] + am2m1*m[i])/(am4m3 + am2m1);
            }
            else {
                t[i] = 0.5 * (m[i-1] + m[i]);
            }
        }

        for (int i=2; i<N-3; i++) {
            double p0 = y[i];
            double p1 = t[i];
            double dx = (x[i+1]-x[i]);
            double p2 = (3*m[i] - 2*t[i] - t[i+1])/dx;
            double p3 = (-2*m[i] + t[i] + t[i+1])/(dx*dx);
            for (int j=0; j<nSubPoints; j++) {
                dx = j*(x[i+1] - x[i])/nSubPoints;
                int iData = (i-2)*nSubPoints+j;
                interpolatedX[iData] = x[i]+dx;
                interpolatedY[iData] = p0 + p1*dx + p2*dx*dx + p3*dx*dx*dx;
            }
        }
    }
    
    public static void main(String[] args) {
        FileReader fileReader;
        try {
            fileReader = new FileReader("B7.dat.scaled");
        }catch(IOException e) {
            throw new RuntimeException("Cannot open B7.dat.scaled, caught IOException: " + e.getMessage());
        }
        ArrayList<Double> xList = new ArrayList<Double>();
        ArrayList<Double> yList = new ArrayList<Double>();
        try {
            BufferedReader bufReader = new BufferedReader(fileReader);
            while (true) {
                String line = bufReader.readLine();
                if (line == null) {
                    break;
                }
                String[] xy = line.split(" +");
                xList.add(Double.parseDouble(xy[0]));
                yList.add(Double.parseDouble(xy[1]));
            }
            fileReader.close();
        } catch(IOException e) {
            throw new RuntimeException("Problem reading d.dat, caught IOException: " + e.getMessage());
        }

        double[] x = new double[xList.size()];
        double[] y = new double[yList.size()];
        for (int i=0; i<x.length; i++) {
            x[i] = xList.get(i);
            y[i] = yList.get(i);
        }
        AkimaSpline fitter = new AkimaSpline(10);
        fitter.init(x.length);
        fitter.doInterpolation(x, y);
        for (int i=0; i<fitter.interpolatedX.length; i++) {
            System.out.println(fitter.interpolatedX[i]+" "+fitter.interpolatedY[i]);
        }
    }

    // hack for main method to prep for doInterpolation
    protected void init(int l) {
        interpolatedX = new double[nSubPoints*(l-1)+1];
        interpolatedY = new double[interpolatedX.length];
    }        

    protected IData processData(IData inputData) {
        double[] originalX = originalDataInfo.getXDataSource().getIndependentData(0).getData();
        double[] originalY = ((DataDoubleArray)inputData).getData();
        doInterpolation(originalX, originalY);
        return interpolatedData;
    }

    protected IEtomicaDataInfo processDataInfo(IEtomicaDataInfo inputDataInfo) {
        originalDataInfo = (DataInfoFunction)inputDataInfo;
        int l = inputDataInfo.getLength();
        if (interpolatedX.length != nSubPoints*(l-1)+1) {
            IEtomicaDataInfo originalXInfo = originalDataInfo.getXDataSource().getIndependentDataInfo(0);
            if (l < 3) throw new RuntimeException("need at least 3 data points");
            if (l != originalXInfo.getLength()) throw new RuntimeException("x and y have different lengths");
            xOutDataInfo = new DataInfoDoubleArray(originalXInfo.getLabel(), originalXInfo.getDimension(), new int[]{nSubPoints*(l-1)+1});
            interpolatedX = new double[nSubPoints*(l-1)+1];
            xOutData = new DataDoubleArray(new int[]{interpolatedX.length}, interpolatedX);
            interpolatedY = new double[interpolatedX.length];
            interpolatedData = new DataFunction(new int[]{interpolatedY.length}, interpolatedY);
            dataInfo = new DataInfoFunction("interpolated "+inputDataInfo.getLabel(), inputDataInfo.getDimension(), this);
        }
        return dataInfo;
    }

    public DataPipe getDataCaster(IEtomicaDataInfo intputDataInfo) {
        return null;
    }
    
    public int getIndependentArrayDimension() {
        return 1;
    }

    public DataDoubleArray getIndependentData(int i) {
        return xOutData;
    }

    public DataInfoDoubleArray getIndependentDataInfo(int i) {
        return xOutDataInfo;
    }

    public DataTag getIndependentTag() {
        return xTag;
    }

    protected DataDoubleArray xOutData;
    protected DataInfoDoubleArray xOutDataInfo;
    protected DataFunction interpolatedData;
    protected DataInfoFunction dataInfo, originalDataInfo;
    protected double[] m, t, x, y, interpolatedX, interpolatedY;
    protected int nSubPoints;
    protected final DataTag xTag;
}
