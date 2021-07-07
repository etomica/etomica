/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.mappedDensity.crystal;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.DataSourceUniform.LimitType;
import etomica.data.histogram.Histogram;
import etomica.data.histogram.HistogramExpanding;
import etomica.data.histogram.HistogramSimple;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.math.DoubleRange;
import etomica.normalmode.CoordinateDefinition;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.space.Vector;
import etomica.units.dimensions.*;

import java.util.Objects;

/**
 * Measured anisotropic singlet density for a crystal using HMA, focusing on variation of a single coordinate variable
 */
public class MeterDensityAnisotropic1CoordHMA implements IDataSource, DataSourceIndependent, AtomLeafAgentManager.AgentSource<Vector> {

    protected final PotentialMaster potentialMaster;
    protected final IteratorDirective id;
    protected DataSourceUniform xDataSourceR;
    protected DataSourceUniform xDataSourceTheta;
    protected DataSourceUniform xDataSourcePhi;

    protected final PotentialCalculationForceSum pc;
    protected final AtomLeafAgentManager<Vector> agentManager;
    protected final Box box;
    protected DataFunction data;
    protected IDataInfo dataInfo;
    protected double Rmax;
    protected final Vector rivector;
    protected final Vector deltaVec;
    protected final Vector rVec;
    protected final Vector RVec;
    protected final Vector riDotB;
    protected final Vector riDotA;
    protected final DataTag tag;
    protected double temperature;
    protected double beta;
    protected CoordinateDefinition latticesite;
    protected double msd;
    protected int nR;
    protected int nTheta;
    protected int nPhi;

    /**
     *
     * @param msd mean-square displacement defining Gaussian reference
     * @param iX index specifying coordinate that is varied. 0 --> r, 1 --> theta, 2 --> phi
     * @param nX number of tabulated values of coordinate specified by iX
     * @param X2 fixed value of first non-iX variable
     * @param X3 fixed value of second non-iX variable
     * @param box
     * @param potentialMaster
     * @param temperature
     * @param latticesite
     */
    public MeterDensityAnisotropic1CoordHMA(double msd, int iX, int nX, double X2, double X3, Box box, PotentialMaster potentialMaster, double temperature, CoordinateDefinition latticesite) {
        this.box = box;
        this.temperature = temperature;
//        this.temperature = 1.0;
        beta = 1./this.temperature;
        this.msd = msd;
        System.out.println("sigma: "+Math.sqrt(msd/3.));
        this.Rmax = Math.sqrt(msd)*2;
        this.latticesite = latticesite;
        this.rivector = box.getSpace().makeVector();
        this.deltaVec = box.getSpace().makeVector();
        this.rVec = box.getSpace().makeVector();
        this.RVec = box.getSpace().makeVector();
        this.riDotB = box.getSpace().makeVector();
        this.riDotA = box.getSpace().makeVector();

        xDataSourceR = new DataSourceUniform("r", Length.DIMENSION);
        xDataSourceR.setTypeMax(LimitType.HALF_STEP);
        xDataSourceR.setTypeMin(LimitType.HALF_STEP);

        xDataSourceTheta = new DataSourceUniform("theta", Angle.DIMENSION);
        xDataSourceTheta.setTypeMax(LimitType.HALF_STEP);
        xDataSourceTheta.setTypeMin(LimitType.HALF_STEP);

        xDataSourcePhi = new DataSourceUniform("theta", Null.DIMENSION);
        xDataSourcePhi.setTypeMax(LimitType.HALF_STEP);
        xDataSourcePhi.setTypeMin(LimitType.HALF_STEP);

        switch(iX) {
            case 0: //r
                xDataSourceR.setXMin(0);
                xDataSourceR.setXMax(Rmax);
                xDataSourceTheta.setXMin(X2);
                xDataSourceTheta.setXMax(X2);
                xDataSourcePhi.setXMin(X3);
                xDataSourcePhi.setXMax(X3);
                nR = nX;
                nTheta = 1;
                nPhi = 1;
                break;
            case 1: //theta
                xDataSourceR.setXMin(X2);
                xDataSourceR.setXMax(X2);
                xDataSourceTheta.setXMin(-1);
                xDataSourceTheta.setXMax(+1);
                xDataSourcePhi.setXMin(X3);
                xDataSourcePhi.setXMax(X3);
                nR = 1;
                nTheta = nX;
                nPhi = 1;
                break;
            case 2: //phi
                xDataSourceR.setXMin(X2);
                xDataSourceR.setXMax(X2);
                xDataSourceTheta.setXMin(X3);
                xDataSourceTheta.setXMax(X3);
                xDataSourcePhi.setXMin(0);
                xDataSourcePhi.setXMax(2*Math.PI);
                nR = 1;
                nTheta = 1;
                nPhi = nX;
                break;
            default: throw new IllegalArgumentException("iX must be 0, 1, or 2");
        }
        xDataSourceR.setNValues(nR);
        xDataSourceTheta.setNValues(nTheta);
        xDataSourcePhi.setNValues(nPhi);

        tag = new DataTag();
        this.potentialMaster = potentialMaster;
        id = new IteratorDirective();
        pc = new PotentialCalculationForceSum();
        agentManager = new AtomLeafAgentManager<Vector>(this, box);
        pc.setAgentManager(agentManager);

        data = new DataFunction(new int[]{nR,nTheta,nPhi});
        dataInfo = new DataInfoFunction("Mapped Average Profile", new CompoundDimension(new Dimension[]{Quantity.DIMENSION, Volume.DIMENSION}, new double[]{1, -1}), this);
        dataInfo.addTag(tag);

    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }
//sim.coordinateDefinition - latticesites

    public final Histogram hl = new HistogramExpanding(1);
    public final Histogram h = new HistogramSimple(200, new DoubleRange(-1e7, 1e7));
    /**
     * Returns the profile for the current configuration.
     */
    public IData getData() {
        pc.reset();
        potentialMaster.calculate(box, id, pc);
        data.E(0);
        double[] y = data.getData();

        IAtomList atoms = box.getLeafList();
        double [][][] ytemporary=new double[xDataSourceR.getNValues()][xDataSourceTheta.getNValues()][xDataSourcePhi.getNValues()];

        double sigma2 = msd / 3;
        double sigma = Math.sqrt(sigma2);
        double q = Math.pow(2 * Math.PI * sigma2, 1.5);

        for (IAtom atom : atoms) {
            rivector.Ev1Mv2(atom.getPosition(), latticesite.getLatticePosition(atom)); // r - R
            box.getBoundary().nearestImage(rivector);
            double ri = Math.sqrt(rivector.squared());//distance of atom from lattice site, r-R
            riDotA.E(Singlet3DmappingDelta0.xyzDotA(rivector, ri, sigma));
            double a1 = -beta * agentManager.getAgent(atom).dot(riDotA);

            double dlnpridr = - ri / sigma2;
            //double rdot = riDotA.dot(rivector) / ri;
            double rdot = Math.sqrt(riDotA.squared());//riDotA is already in the direction of rivector
            double a2 = rdot * dlnpridr;

            for (int i = 0; i < nR; i++) {
                double r = ((DataDoubleArray) xDataSourceR.getData()).getData()[i];
                double pr = (Math.exp(-r * r / (2 * sigma2)));

                for (int j = 0; j < nTheta; j++) {
                    double costheta = ((DataDoubleArray) xDataSourceTheta.getData()).getData()[j];
                    double sintheta = Math.sqrt(1.0 - costheta*costheta);

                    for (int k = 0; k < nPhi; k++) {
                        double phi = ((DataDoubleArray) xDataSourcePhi.getData()).getData()[k];

                        //point where dens measured, relative to lat site;  delta - R
                        deltaVec.setX(0, r * sintheta * Math.cos(phi));
                        deltaVec.setX(1, r * sintheta * Math.sin(phi));
                        deltaVec.setX(2, r * costheta);

                        rVec.Ev1Mv2(rivector, deltaVec); // r - delta
                        riDotB.E(Singlet3DmappingDelta0.xyzDotB(rVec, Math.sqrt(rVec.squared())));

                        double a3 = -beta * agentManager.getAgent(atom).dot(riDotB);
//                        double r2 = rVec.squared()/sigma2;
                        ytemporary[i][j][k] += (pr / q) + (a1 + a2) * pr + a3;
                        if(i==0 && j==0 && k==0 && Objects.equals(atom, atoms.get(0))) {
//                            System.out.println((a1 + a2) * pr + a3+", "+pr+", "+a1*pr+", "+a2*pr+", "+a3);
                        }
//                        if  (k==0) {
//                            h.addValue( - agentManager.getAgent(atom).dot(vel) + (temperature * rdot * dlnpridr));
//                            hl.addValue(Math.log(Math.abs(- agentManager.getAgent(atom).dot(vel) + (temperature * rdot * dlnpridr))));
//                        }
                        if (Double.isNaN(ytemporary[i][j][k])) {
                            Singlet3DmappingDelta0.xyzDot(rivector, deltaVec, sigma);
                            System.out.println("oops");
                        }
                    }
                }
            }
            //debugging
            double myr = Math.sqrt(rVec.squared()); // |ri - delta|
            //if(i==0 && j==0 && k==0 && Math.abs((a1+a2)*pr+a3)>1e6) System.out.println((ri/sigma)+", "+myr+", "+a1*pr+", "+a2*pr+", "+a3+", "+((a1+a2)*pr+a3)+", "+ Math.sqrt(riDotA.squared())+", "+Math.sqrt(riDotB.squared()));
//                        rVec.E(rivector); // atom position relative to lattice site
//                        rVec.normalize(); //unit vector from lattice site to atom
//                        double Fr = agentManager.getAgent(atom).dot(rVec);//magnitude of force in direction from lattice site to atom
//                        rVec.TE(Fr);//force in direction of r-R
//                        rVec.ME(agentManager.getAgent(atom));//force in direction perpendicular to r-R
//                        double Fperp = Math.sqrt(rVec.squared());
//                        //if(new java.util.Random().nextDouble() < 0.01) System.out.println(Fperp+", "+Fr+", "+temperature*dlnpridr+", "+(Fr-temperature*dlnpridr)+", "+(Fperp/Fr));

            //end debugging

        }//loop over atoms
        int n = 0;
        for (int i = 0; i < nR; i++) {
            for (int j = 0; j < nTheta; j++) {
                for (int k = 0; k < nPhi; k++) {
                    y[n] = ytemporary[i][j][k] / box.getLeafList().size();
                    n = n + 1;
                }
            }
        }
         return data;
    }

    public DataDoubleArray getIndependentData(int i) {
        if(i==0){return (DataDoubleArray) xDataSourceR.getData();}
        if(i==1){return (DataDoubleArray) xDataSourceTheta.getData();}
        if(i==2){return (DataDoubleArray) xDataSourcePhi.getData();}
        throw new RuntimeException();
    }
    public DataInfoDoubleArray getIndependentDataInfo(int i) {

        if(i==0){return (DataInfoDoubleArray) xDataSourceR.getDataInfo();}
        if(i==1){return (DataInfoDoubleArray) xDataSourceTheta.getDataInfo();}
        if(i==2){return (DataInfoDoubleArray) xDataSourcePhi.getDataInfo();}
        throw new RuntimeException();
    }
    public int getIndependentArrayDimension() {
        return 3;
    }

    public DataTag getIndependentTag(){throw new RuntimeException();}
    public DataTag getIndependentTag(int i) {

        if(i==0){return xDataSourceR.getTag();}
        if(i==1){return xDataSourceTheta.getTag();}
        if(i==2){return xDataSourcePhi.getTag();}
        throw new RuntimeException();
    }
    /**
     * @return Returns the box.
     */
    public Box getBox() {
        return box;
    }

    public DataSourceUniform getXDataSource(int i) {
        if(i==0)return xDataSourceR;
        if(i==1)return xDataSourceTheta;
        if(i==2)return xDataSourcePhi;
        throw new RuntimeException();

    }

    @Override
    public Vector makeAgent(IAtom a, Box agentBox) {
        return box.getSpace().makeVector();
    }

    @Override
    public void releaseAgent(Vector agent, IAtom atom, Box agentBox) {

    }
}
