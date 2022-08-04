package etomica.mappedRdf;

import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.potential.compute.PotentialCompute;
import etomica.units.dimensions.Null;

/**
 * Calculates pair distribution using mapped averaging
 */
public class MeterMappedRdf implements IDataSource, DataSourceIndependent {

    protected Box box;
    protected DataFunction data;
    private IDataInfo dataInfo;
    protected DataDoubleArray rData;
    protected final DataSourceUniform xDataSource;
    protected double xMax;
    protected final DataTag tag;
    protected final PotentialCompute potentialMaster;
    protected final PotentialCallbackMappedRdf pc;
    protected double density;
    protected double rcforHandfinmap;

    public MeterMappedRdf(double rcforHandfinmap, PotentialCompute potentialMaster, Box box, int nbins, double density) {
        this.box = box;
        this.density = density;
        this.potentialMaster = potentialMaster;
        this.rcforHandfinmap = rcforHandfinmap;

        pc = new PotentialCallbackMappedRdf(rcforHandfinmap, box, nbins, potentialMaster);

        xDataSource = pc.getXDataSource();

        rData = (DataDoubleArray) xDataSource.getData();
        data = new DataFunction(new int[]{rData.getLength()});
        dataInfo = new DataFunction.DataInfoFunction("g(r)", Null.DIMENSION, this);

        tag = new DataTag();
        dataInfo.addTag(tag);
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    /**
     * Zero's out the RDF sum tracked by this meter.
     */
    public void reset() {
        rData = (DataDoubleArray) xDataSource.getData();
        xMax = xDataSource.getXMax();
        pc.getXDataSource().setNValues(rData.getLength());
        pc.getXDataSource().setXMax(xMax);
        data = new DataFunction(new int[]{rData.getLength()});
        dataInfo = new DataFunction.DataInfoFunction("g(r)", Null.DIMENSION, this);
        dataInfo.addTag(tag);
    }

    /**
     * Returns the RDF, averaged over the calls to actionPerformed since the
     * meter was reset or had some parameter changed (xMax or # of bins).
     */
    public IData getData() {
        if (rData != xDataSource.getData() ||
                data.getLength() != rData.getLength() ||
                xDataSource.getXMax() != xMax) {
            reset();
        }

        // just to get the forces
        potentialMaster.computeAll(true);

        // now do our mapped rdf
        pc.reset();
        potentialMaster.computeAll(false, pc);
        potentialMaster.computeAll(false, pc);

        final double[] y = data.getData();

        double[] r = rData.getData();
        double[] gSum = pc.getGSum();
        double vol = box.getBoundary().volume();
        int numAtoms = box.getLeafList().size();

        for (int i = 0; i < r.length; i++) {
//            double vShell = space.sphereVolume(r[i]+dx2)-space.sphereVolume(r[i]-dx2);
            // y[i] = (gR[i]*callCount+gSum[i])*01/(norm*vShell);

            y[i] = numAtoms * (numAtoms - 1) / (vol * vol) + (gSum[i]);
            //    y[i] =  (gSum[i] / (norm)) ;

        }
        //       System.out.println(y[10]);
        return data;
    }

    public DataSourceUniform getXDataSource() {
        return xDataSource;
    }

    public DataDoubleArray getIndependentData(int i) {
        return (DataDoubleArray) xDataSource.getData();
    }

    public DataDoubleArray.DataInfoDoubleArray getIndependentDataInfo(int i) {
        return (DataDoubleArray.DataInfoDoubleArray) xDataSource.getDataInfo();
    }

    public DataTag getIndependentTag() {
        return xDataSource.getTag();
    }

    public int getIndependentArrayDimension() {
        return 1;
    }

    /**
     * @return Returns the box.
     */
    public Box getBox() {
        return box;
    }

    /**
     * @param box The box to set.
     */
    public void setBox(Box box) {
        pc.setBox(box);
        this.box = box;
    }

    public PotentialCallbackMappedRdf getPotentialCallback() {
        return pc;
    }

}
