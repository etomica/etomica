package etomica.spin.heisenberg;

import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.units.dimensions.Null;


public class MeterMappedAveragingCV implements IDataSource {
    protected final DataDoubleArray data;
    protected final DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected double bJ;
    protected double bt;
    protected final int L;
    protected final int N;
    protected final int[][] nbrs = new int[3][4];
    protected final double[] sintiMtj1;
    protected final double[] sintiMtj2;
    protected final double[] sintiMtj3;

    protected final double[] costiMtj1;
    protected final double[] sin2tiMtj1;
    protected final double[] sin2tiMtjtk1;

    protected double temperature;
    private Box box;
    protected final Space space;
    protected final PotentialMaster potentialMaster;

    public MeterMappedAveragingCV(Box box, double temperature, double interactionS, PotentialMaster potentialMaster) {
        this.box = box;
        this.space = box.getSpace();
        this.temperature = temperature;
        this.potentialMaster = potentialMaster;
        bt = 1 / temperature;
        bJ = interactionS * bt;

        data = new DataDoubleArray(2);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{2});
        tag = new DataTag();
        dataInfo.addTag(tag);
        N = box.getLeafList().size();
        L = (int) Math.round(Math.sqrt(N));
        sintiMtj1 = new double[N];
        sintiMtj2 = new double[N];
        sintiMtj3 = new double[N];
        costiMtj1 = new double[N];
        sin2tiMtj1 = new double[N];
        sin2tiMtjtk1 = new double[N];

    }


    public IData getData() {
        double[] x = data.getData();
        IAtomList leafList = box.getLeafList();
        //double term1 = 0, term2 = 0, term2a = 0, term2b = 0, term2c = 0, term2d = 0, term3a = 0, term3b = 0, term4a = 0, term4b= 0, term5 = 0;

        //compute neighbor sums
        for (int i = 0; i < N; i++) {
            double costi = ((IAtomOriented) leafList.get(i)).getOrientation().getDirection().getX(0);
            double sinti = ((IAtomOriented) leafList.get(i)).getOrientation().getDirection().getX(1);

            sintiMtj1[i] = 0;
            costiMtj1[i] = 0;
            sin2tiMtj1[i] = 0;
            sin2tiMtjtk1[i] = 0;
            sintiMtj2[i] = 0;
            sintiMtj3[i] = 0;

            getNeighbors(i);

            //1st-neighbor loop
            for (int j = 0; j < 4; j++) {
                double costj = ((IAtomOriented) leafList.get(nbrs[0][j])).getOrientation().getDirection().getX(0);
                double sintj = ((IAtomOriented) leafList.get(nbrs[0][j])).getOrientation().getDirection().getX(1);

                double s = sinti * costj - costi * sintj;//sin(ti-tj)
                double c = costi * costj + sinti * sintj;//cos(ti-tj)
                sintiMtj1[i] += s;
                costiMtj1[i] += c;
                sin2tiMtj1[i] += 2 * s * c;//sin(2(ti-tj))

                for (int k = 0; k < j; k++) {
                    double costk = ((IAtomOriented) leafList.get(nbrs[0][k])).getOrientation().getDirection().getX(0);
                    double sintk = ((IAtomOriented) leafList.get(nbrs[0][k])).getOrientation().getDirection().getX(1);

//                    sin(2ti-tj-tk) = sin(2ti)Cos(tj+tk) - cos(2ti)*sin(tj+tk)
                    double costjPtk = costj * costk - sintj * sintk;
                    double sintjPtk = sintj * costk + costj * sintk;
                    sin2tiMtjtk1[i] += 2 * sinti * costi * costjPtk - (2 * costi * costi - 1) * sintjPtk;
                }
            }

            //2nd-neighbor loop
            for (int j2 = 0; j2 < 4; j2++) {
                double costj2 = ((IAtomOriented) leafList.get(nbrs[1][j2])).getOrientation().getDirection().getX(0);
                double sintj2 = ((IAtomOriented) leafList.get(nbrs[1][j2])).getOrientation().getDirection().getX(1);
                sintiMtj2[i] += sinti * costj2 - costi * sintj2;
            }

            //3rd-neighbor loop
            for (int j3 = 0; j3 < 4; j3++) {
                double costj3 = ((IAtomOriented) leafList.get(nbrs[2][j3])).getOrientation().getDirection().getX(0);
                double sintj3 = ((IAtomOriented) leafList.get(nbrs[2][j3])).getOrientation().getDirection().getX(1);
                sintiMtj3[i] += sinti * costj3 - costi * sintj3;
            }
        }

        double sumOverI = 0;
        double sumVar = 0;
        for (int i = 0; i < N; i++) {
            double costi = ((IAtomOriented) leafList.get(i)).getOrientation().getDirection().getX(0);
            double sinti = ((IAtomOriented) leafList.get(i)).getOrientation().getDirection().getX(1);

            getNeighbors(i);

            double sum = 0;
            for (int j = 0; j < 4; j++) {
                double costj = ((IAtomOriented) leafList.get(nbrs[0][j])).getOrientation().getDirection().getX(0);
                double sintj = ((IAtomOriented) leafList.get(nbrs[0][j])).getOrientation().getDirection().getX(1);
                sum += (costi * costj + sinti * sintj) * sintiMtj1[nbrs[0][j]];
            }

            sum += -costiMtj1[i] * sintiMtj1[i];
            sum *= bJ * bJ * bJ * 0.25 * sintiMtj1[i]; //multiply by -(bJ^3)/4 (fi/J)
            sum *= 2; //phi contribution is the same, so add it in by doubling

            //bThetadot^beta_beta
            double tbb = 0.125 * (sin2tiMtj1[i] + 2 * sin2tiMtjtk1[i] - 4 * sintiMtj2[i] - 2 * sintiMtj3[i]);

            sum += -bJ * bJ * bJ * sintiMtj1[i] * tbb;

            sumOverI += sum;

            sumVar += sintiMtj1[i]*sintiMtj1[i];

        }

        x[0] = N*bJ*bJ + sumOverI + 0.25*bJ*bJ*bJ*bJ*sumVar*sumVar;
        x[1] = -0.5*bJ*bJ*sumVar;//part of variance contribution. The average of this is squared and subtracted from the x[0] average

        return data;
    }


    public DataTag getTag() {
        return tag;
    }


    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    //creates arrays of first, second, and third neighbors from given site
    private void getNeighbors(int i) {
        boolean right = i % L == L - 1;
        boolean left = i % L == 0;
        boolean top = i < L;
        boolean bottom = i >= (N - L);
        boolean right2 = i % L >= L - 2;
        boolean left2 = i % L <= 1;
        boolean top2 = i < 2 * L;
        boolean bottom2 = i >= N - 2 * L;

        nbrs[0][0] = i + 1;
        nbrs[0][1] = i - 1;
        nbrs[0][2] = i - L;
        nbrs[0][3] = i + L;
        nbrs[1][0] = i + L + 1;//right bottom
        nbrs[1][1] = i + L - 1;//left bottom
        nbrs[1][2] = i - L - 1;//left top
        nbrs[1][3] = i - L + 1;//right top

        nbrs[2][0] = i + 2;
        nbrs[2][1] = i - 2;
        nbrs[2][2] = i - 2 * L;
        nbrs[2][3] = i + 2 * L;

        if (right) {
            nbrs[0][0] -= L;
            nbrs[1][0] -= L;
            nbrs[1][3] -= L;
        } else if (left) {
            nbrs[0][1] += L;
            nbrs[1][1] += L;
            nbrs[1][2] += L;
        }
        if (top) {
            nbrs[0][2] += N;
            nbrs[1][2] += N;
            nbrs[1][3] += N;
        } else if (bottom) {
            nbrs[0][3] -= N;
            nbrs[1][0] -= N;
            nbrs[1][1] -= N;
        }

        if (right2) {
            nbrs[2][0] -= L;
        } else if (left2) {
            nbrs[2][1] += L;
        }
        if (top2) {
            nbrs[2][2] += N;
        } else if (bottom2) {
            nbrs[2][3] -= N;
        }


    }


}
