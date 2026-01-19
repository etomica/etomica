package etomica.data.histogram;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.math.DoubleRange;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;

import java.util.List;
import java.util.Set;

public class HistogramVectorSource {
    protected double deltaX, deltaY, deltaZ;
    protected Box box;
    protected double sum;
    protected long[][][] counts;
    protected double[][][] energyHist;
    protected double[][][] histogram;
    protected double[] xValues, yValues, zValues;
    protected double xMin, yMin, zMin;
    protected double xMax, yMax, zMax, LDnBin, L;
    protected ISpecies species;
    protected Set<List<Integer>> setUnique;
    protected Integer[] set = new Integer[3];
    protected int nBinsX, nBinsY, nBinsZ;
    protected long numsteps;
    protected final boolean rectangular;


    public HistogramVectorSource(int nX, int nY, int nZ, DoubleRange rangeX, DoubleRange rangeY, DoubleRange rangeZ, ISpecies species, Box box, long numsteps){
        nBinsX = nX;
        nBinsY = nY;
        nBinsZ = nZ;
        this.rectangular = box.getBoundary().isRectangular();
        counts = new long[nX][nY][nZ];
        // histogram = new double[nX][nY][nZ];
        energyHist = new double[nX][nY][nZ];
        xValues = new double[nX];
        yValues = new double[nY];
        zValues = new double[nZ];
        setXRange(rangeX);
        setYRange(rangeY);
        setZRange(rangeZ);
        setBox(box);
        setSpecies(species);
        setL(box);
        setNumSteps(numsteps);
        //numsteps = box.ge
    }

    private void setBox(Box box) {
        this.box = box;
    }

    private void setNumSteps(long steps){this.numsteps = steps;}

    //set box length
    private void setL(Box box){
        L = box.getBoundary().getBoxSize().getX(0);
    }

    public Box getBox(){return box;}

    //set range
    public void setRange(DoubleRange rangeX, DoubleRange rangeY, DoubleRange rangeZ) {
        xMin = rangeX.minimum();
        xMax = rangeX.maximum();
        yMin = rangeY.minimum();
        yMax = rangeY.maximum();
        zMin = rangeZ.minimum();
        zMax = rangeZ.maximum();
        deltaX = (xMax-xMin)/nBinsX;
        deltaY = (yMax-yMin)/nBinsY;
        deltaZ = (zMax-zMin)/nBinsZ;
        reset();
    }

    public void reset() {
        sum = 0;
        for (int i = 0; i < nBinsX; i++) {
            for (int j = 0; j < nBinsY; j++) {
                for (int k = 0; k < nBinsZ; k++) {
                    counts[i][j][k] = 0;
                    xValues[i] = xMin + deltaX * (i + 0.5);
                    yValues[j] = yMin + deltaY * (j + 0.5);
                    zValues[k] = zMin + deltaZ * (k + 0.5);
                }
            }
        }

    }

    // add value into the histogram
   /* public void addValue(double x, double y, double z, double energy) {
        //   if( x > -3 && y > -3 && z > -3){
        int i = (int) ((x + L / 2) / (deltaX));
        int j = (int) ((y + L / 2) / (deltaY));
        int k = (int) ((z + L / 2) / (deltaZ));

        // Ensure indices are within bounds
        if(i > nBinsX - 1 && i<0  ){
            throw  new RuntimeException("error in i value");
        }

       /* i = Math.min(i, nBinsX - 1);
        j = Math.min(j, nBinsX - 1);
        k = Math.min(k, nBinsX - 1);

        // Increment the count for the appropriate cube
        counts[i][j][k]++;
        energyHist[i][j][k] += energy;
        // }
    }*/

    public void addValue(double x, double y, double z, double energy){

        Vector s = box.getSpace().makeVector();
        Vector r = new Vector3D(x, y, z);
        s.E(r);

        if (box.getBoundary().isRectangular()) {
            Vector bs = box.getBoundary().getBoxSize();
            s.DE(bs);
        } else {
            Tensor hInv = box.getBoundary().getHInv();
            hInv.transform(s);
        }

        double u = s.getX(0)+ 0.5;
        double v = s.getX(1)+ 0.5;
        double w = s.getX(2)+ 0.5;

        int i = ((int)(u * nBinsX) % nBinsX + nBinsX ) % nBinsX;
        int j = (int)(v * nBinsY) % nBinsY;
        int k = (int)(w * nBinsZ) % nBinsZ;

        counts[i][j][k]++;
        energyHist[i][j][k] += energy;
    }

    //set range of histogram
    public void setXRange(DoubleRange xRange) {
        xMin = xRange.minimum();
        xMax = xRange.maximum();
        deltaX = (xMax - xMin) / nBinsX;
        reset();
    }

    public void setYRange(DoubleRange yRange) {
        yMin = yRange.minimum();
        yMax = yRange.maximum();
        deltaY = (yMax - yMin) / nBinsY;
        reset();
    }

    public void setZRange(DoubleRange zRange) {
        zMin = zRange.minimum();
        zMax = zRange.maximum();
        deltaZ = (zMax - zMin) / nBinsZ;
        reset();
    }

    public double[][][] getEnergyHist() {
       /* for (int i = 0; i < nBinsX; i++) {
            for (int j = 0; j < nBinsY; j++) {
                for (int k = 0; k < nBinsZ; k++) {
                    if (counts[i][j][k]  != 0){
                        System.out.println(energyHist[i][j][k]   +" "+counts[i][j][k]);
                        sum = energyHist[i][j][k]   / counts[i][j][k];
                        histogram[i][j][k] = sum;
                    }else {
                        histogram[i][j][k] = 0;
                    }

                }
            }
        }*/
        return energyHist;
    }

    public void setSpecies(ISpecies species){
        this.species = species;
    }

    public ISpecies getSpecies(){
        return species;
    }
    public long getSteps(){return numsteps;}

    public long[][][] getCounts(){return counts;}



}
