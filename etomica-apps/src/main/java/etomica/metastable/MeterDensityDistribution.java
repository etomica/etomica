package etomica.metastable;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.data.AccumulatorHistogram;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.space.Space;
import etomica.units.dimensions.CompoundDimension;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Quantity;
import etomica.units.dimensions.Volume;
import etomica.data.histogram.HistogramDiscrete;

public class MeterDensityDistribution implements IEtomicaDataSource {

    protected AccumulatorHistogram histogram; 
    protected final Box box;
    protected final int[][][] counts;
    protected final Vector p2;
    protected final Vector subBox, shift;
    protected final DataDouble subData;
    
    public MeterDensityDistribution(Space space, Box box, int nSubBoxes) {
        this.box = box;
        histogram = new AccumulatorHistogram(new HistogramDiscrete(1e-10));
        Dimension densityDim = new CompoundDimension(new Dimension[]{Quantity.DIMENSION,Volume.DIMENSION},new double[]{1,-1});
        histogram.putDataInfo(new DataInfoDouble("density", densityDim));
        subData = new DataDouble();
        counts = new int[nSubBoxes][nSubBoxes][nSubBoxes];
        p2 = space.makeVector();
        subBox = space.makeVector();
        shift = space.makeVector();
    }
    
    public void reset() {
        histogram.reset();
    }

    public IData getData() {
        IAtomList atoms = box.getLeafList();
        Vector L = box.getBoundary().getBoxSize();
        int nSubBoxes = counts.length;
        subBox.Ea1Tv1(1.0/nSubBoxes, L);
        double vSubBox = box.getBoundary().volume()/Math.pow(nSubBoxes, 3);
        shift.E(0);
        for (int is=0; is<2; is++) {
            shift.setX(0, 0.5*is*subBox.getX(0));
            for (int js=0; js<2; js++) {
                shift.setX(1, 0.5*js*subBox.getX(1));
                for (int ks=0; ks<2; ks++) {
                    shift.setX(2, 0.5*ks*subBox.getX(2));
                    
                    for (int i=0; i<nSubBoxes; i++) {
                        for (int j=0; j<nSubBoxes; j++) {
                            for (int k=0; k<nSubBoxes; k++) {
                                counts[i][j][k] = 0;
                            }
                        }
                    }
                    for (int i=0; i<atoms.getAtomCount(); i++) {
                        p2.E(atoms.getAtom(i).getPosition());
                        p2.PE(shift);
                        p2.PEa1Tv1(0.5, L);
                        p2.DE(subBox);
                        int idx0 = (int)p2.getX(0);
                        idx0 = idx0 < nSubBoxes ? idx0 : 0;
                        int idx1 = (int)p2.getX(1);
                        idx1 = idx1 < nSubBoxes ? idx1 : 0;
                        int idx2 = (int)p2.getX(2);
                        idx2 = idx2 < nSubBoxes ? idx2 : 0;
                        counts[idx0][idx1][idx2]++;
                    }
                    for (int i=0; i<nSubBoxes; i++) {
                        for (int j=0; j<nSubBoxes; j++) {
                            for (int k=0; k<nSubBoxes; k++) {
                                subData.x = counts[i][j][k] / vSubBox;
                                histogram.putData(subData);
                            }
                        }
                    }
                }
            }
        }
        return histogram.getData();
    }

    public DataTag getTag() {
        return histogram.getTag();
    }

    public IEtomicaDataInfo getDataInfo() {
        return histogram.getDataInfo();
    }

}
