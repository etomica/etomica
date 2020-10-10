package etomica.mappedRdf;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomPair;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.box.storage.Tokens;
import etomica.box.storage.VectorStorage;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

/**
 * Calculates pair distribution using mapped averaging
 */
public class MeterMappedRdf implements IDataSource, DataSourceIndependent, AtomLeafAgentManager.AgentSource<Vector> {

    protected final PotentialCalculationForceSum pcForce;
    protected double density;
    protected double rcforHandfinmap;

    public MeterMappedRdf(double rcforHandfinmap, Space space, PotentialMaster potentialMaster, Box box, int nbins, double density) {
        this.space = space;
        this.box = box;
        this.density = density;
        this.potentialMaster = potentialMaster;
        this.rcforHandfinmap = rcforHandfinmap;

        VectorStorage forces = box.getAtomStorage(Tokens.vectorsDefault());
        pcForce = new PotentialCalculationForceSum(forces);
        pc = new PotentialCalculationMappedRdf(rcforHandfinmap, space, box, nbins, forces);

        xDataSource = pc.getXDataSource();

        rData = (DataDoubleArray) xDataSource.getData();
        data = new DataFunction(new int[]{rData.getLength()});
        dataInfo = new DataFunction.DataInfoFunction("g(r)", Null.DIMENSION, this);

        allAtoms = new IteratorDirective();

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
        pcForce.reset();
        pc.reset();
        potentialMaster.calculate(box, allAtoms, pcForce);
        long numAtoms = box.getLeafList().size();

        potentialMaster.calculate(box, allAtoms, pc);
        AtomPair foo = new AtomPair();
        for (int i = 0; i < numAtoms; i++) {

            foo.atom0 = box.getLeafList().get(i);
            for (int j = i + 1; j < numAtoms; j++) {
                foo.atom1 = box.getLeafList().get(j);
                pc.doCalculation(foo, null);
            }
        }

        final double[] y = data.getData();

        double[] r = rData.getData();
        double[] gSum = pc.getGSum();
        double vol = box.getBoundary().volume();
        //      System.out.println("metervol " + box.getBoundary().volume());

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

    public String getName() {
        return name;
    }

    public PotentialCalculationMappedRdf getPotentialCalculation() {
        return pc;
    }

    public void setName(String name) {
        this.name = name;
    }


    protected Box box;
    protected final Space space;
    protected DataFunction data;
    private IDataInfo dataInfo;
    protected DataDoubleArray rData;
    protected IteratorDirective allAtoms;
    protected final DataSourceUniform xDataSource;
    protected double xMax;
    private String name;
    protected final DataTag tag;
    protected final PotentialMaster potentialMaster;
    protected final PotentialCalculationMappedRdf pc;

    public Vector makeAgent(IAtom a, Box agentBox) {
        return space.makeVector();
    }

    public void releaseAgent(Vector agent, IAtom atom, Box agentBox) {
    }

}
