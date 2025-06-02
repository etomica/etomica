package etomica.data.histogram;

import etomica.box.Box;
import etomica.math.DoubleRange;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.space.Vector;
import etomica.species.ISpecies;

import java.awt.event.ActionEvent;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class HistogramVectorSimple {
    protected double deltaX, deltaY, deltaZ;
    protected Box box;
    protected long sum;
    protected long[][][] counts;
    protected double[][][] histogram;
    protected double[] xValues, yValues, zValues;
    protected double xMin, yMin, zMin;
    protected double xMax, yMax, zMax, LDnBin, L;
    protected ISpecies species;
    protected Set<List<Integer>> setUnique;
    protected Integer[] set = new Integer[3];
    protected int nBinsX, nBinsY, nBinsZ;
    protected int numsteps;


    public HistogramVectorSimple(int nX, int nY, int nZ, DoubleRange rangeX, DoubleRange rangeY, DoubleRange rangeZ, ISpecies species, Box box, int numsteps){
        nBinsX = nX;
        nBinsY = nY;
        nBinsZ = nZ;
        counts = new long[nX][nY][nZ];
        histogram = new  double[nX][nY][nZ];
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

    private void setNumSteps(int steps){this.numsteps = steps;}

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

   /* public void addValue(double x, double y, double z) {
        if (x >= xMin && x <= xMax && y >= yMin && y <= yMax && z >= zMin && z <= zMax) {
            int i= (int) Math.floor(x/deltaX);
            int j= (int) Math.floor(y/deltaY);
            int k= (int) Math.floor(z/deltaZ);
         /*   int i = (int) Math.ceil((x - xMin) / deltaX);
            int j = (int) Math.ceil((y - yMin) / deltaY);
            int k = (int) Math.ceil((z - zMin) / deltaZ);
            if(x<0 ){
                i = (int) ((i-xMin)/deltaX)-1;
            } else {
                i = (int) ((i-xMin)/deltaX)+1;
            }
            if(y<0){
                j = (int) ((j-yMin)/deltaY)-1;
            }else {
                j = (int) ((j-yMin)/deltaY)+1;
            }
            if(z<0){
                k = (int) ((k-zMin)/deltaZ)-1;
            }else {
                k = (int) ((k-zMin)/deltaZ)+1;
            }
            i = (int) ((i-xMin)/deltaX);
            j = (int) ((j-yMin)/deltaY);
            k = (int) ((k-zMin)/deltaZ);
            if (i == nBinsX) i--;
            if (j == nBinsY) j--;
            if (k == nBinsZ) k--;
            counts[i][j][k]++;
            sum++;
        }
    }*/

    // add value into the histogram
    public void addValue(double x, double y, double z) {
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
        k = Math.min(k, nBinsX - 1);*/

            // Increment the count for the appropriate cube
            counts[i][j][k]++;
       // }

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

    public double[][][] getHistogram() {
            for (int i = 0; i < nBinsX; i++) {
                for (int j = 0; j < nBinsY; j++) {
                    for (int k = 0; k < nBinsZ; k++) {
                        histogram[i][j][k] = counts[i][j][k] / (  deltaX * deltaY * deltaZ );
                    }
                }
            }
        return histogram;
    }

    public void setSpecies(ISpecies species){
        this.species = species;
    }

    public ISpecies getSpecies(){
        return species;
    }
    public int getSteps(){return numsteps;}

    public long[][][] getCounts(){return counts;}

    // molecule coordinates added during simulation
    public void actionPerformed(){
        ISpecies species = getSpecies();
        box = getBox();
        IMoleculeList molecules = box.getMoleculeList(species);
        for(int i=0; i<molecules.size(); i++){
            IMolecule molecule = molecules.get(i);
            //for(int j=0; j<molecule.getChildList().size(); j++){
                Vector vec = molecule.getChildList().get(0).getPosition();
           // System.out.println(vec);
            //System.exit(1);
                addValue(vec.getX(0),vec.getX(1), vec.getX(2) );
          //  }
        }
    }

}
