package etomica.spin.heisenberg_interacting.heisenberg;

import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
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

    public MeterMappedAveragingCV(Simulation sim, double temperature, double interactionS, PotentialMaster potentialMaster) {
        this.box = sim.getBox(0);
        this.space = sim.getSpace();
        this.temperature = temperature;
        this.potentialMaster = potentialMaster;
        bt = 1 / temperature;
        bJ = interactionS * bt;

        data = new DataDoubleArray(2);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{2});
        tag = new DataTag();
        dataInfo.addTag(tag);
        N = box.getLeafList().getAtomCount();
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
        int nM = leafList.getAtomCount();




        for (int i = 0; i < nM; i++) {
            double costi = ((IAtomOriented) leafList.getAtom(i)).getOrientation().getDirection().getX(0);
            double sinti = ((IAtomOriented) leafList.getAtom(i)).getOrientation().getDirection().getX(1);

            sintiMtj1[i] = 0;
            costiMtj1[i] = 0;
            sin2tiMtj1[i] = 0;
            sin2tiMtjtk1[i] = 0;
            getNeighbors(i);
            for (int j = 0; j < 4; j++) {
                double costj = ((IAtomOriented) leafList.getAtom(nbrs[0][j])).getOrientation().getDirection().getX(0);
                double sintj = ((IAtomOriented) leafList.getAtom(nbrs[0][j])).getOrientation().getDirection().getX(1);

                double s = sinti * costj - costi * sintj;
                double c = costi * costj + sinti * sintj;
                sintiMtj1[i] += s;
                costiMtj1[i] += c;
                sin2tiMtj1[i] += 2 * s * c;

                for (int k = 0; k < j; k++) {
                    double costk = ((IAtomOriented) leafList.getAtom(nbrs[0][k])).getOrientation().getDirection().getX(0);
                    double sintk = ((IAtomOriented) leafList.getAtom(nbrs[0][k])).getOrientation().getDirection().getX(1);

//                    sin(2ti-tj-tk) = sin(2ti)Cos(tj+tk) - cos(2ti)*sin(tj+tk)
//
                    double costjPtk = costj * costk - sintj * sintk;
                    double sintjPtk = sintj * costk + costj * sintk;
                    sin2tiMtjtk1[i] += 2 * sinti * costi * costjPtk - (2 * costi * costi - 1) * sintjPtk;


                }


            }

            for (int j2 = 0; j2 < 4; j2++) {
                double costj2 = ((IAtomOriented) leafList.getAtom(nbrs[1][j2])).getOrientation().getDirection().getX(0);
                double sintj2 = ((IAtomOriented) leafList.getAtom(nbrs[1][j2])).getOrientation().getDirection().getX(1);
                sintiMtj2[i] += sinti * costj2 - costi * sintj2;
            }


            for (int j3 = 0; j3 < 4; j3++) {
                double costj3 = ((IAtomOriented) leafList.getAtom(nbrs[2][j3])).getOrientation().getDirection().getX(0);
                double sintj3 = ((IAtomOriented) leafList.getAtom(nbrs[2][j3])).getOrientation().getDirection().getX(1);
                sintiMtj3[i] += sinti * costj3 - costi * sintj3;
            }

        }
        double sumOverI = 0;
        double sumVar = 0;
        for (int i = 0; i < N; i++) {
            double costi = ((IAtomOriented) leafList.getAtom(i)).getOrientation().getDirection().getX(0);
            double sinti = ((IAtomOriented) leafList.getAtom(i)).getOrientation().getDirection().getX(1);


            double sum = 0;
            getNeighbors(i);
            for (int j = 0; j < 4; j++) {
                double costj = ((IAtomOriented) leafList.getAtom(nbrs[0][j])).getOrientation().getDirection().getX(0);
                double sintj = ((IAtomOriented) leafList.getAtom(nbrs[0][j])).getOrientation().getDirection().getX(1);
                sum += (costi * costj + sinti * sintj) * sintiMtj1[nbrs[0][j]];

                sum += sintiMtj1[nbrs[0][j]] * (costi * costj + sinti * sintj);
            }
//            sum += costiMtj1[i] * sintiMtj1[i];//phii
//            sum -= sintiMtj1[i] * costiMtj1[i]; // come from a different line, cancel with phii
            sum *= bJ * bJ * bJ * 0.25 * sintiMtj1[i];


            double tbb = 0.125 * (sin2tiMtj1[i] + 2 * sin2tiMtjtk1[i] - 4 * sintiMtj2[i] - 2 * sintiMtj3[i]);
            sum -= bJ * bJ * bJ * sintiMtj1[i] * tbb;
            sum += bJ * bJ * (1 - sintiMtj1[i] * sintiMtj1[i]);

            sum += bJ*bJ*sintiMtj1[i]*sintiMtj1[i];
            sumOverI += sum;

            sumVar += sintiMtj1[i]*sintiMtj1[i];

        }

        x[0] = sumOverI-0.25*bJ*bJ*sumVar*sumVar;
        x[1] = 0.5*bJ*sumVar;

        return data;
    }


    public DataTag getTag() {
        return tag;
    }


    public IDataInfo getDataInfo() {
        return dataInfo;
    }

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
