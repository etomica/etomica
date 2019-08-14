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
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.math.SpecialFunctions;
import etomica.normalmode.CoordinateDefinition;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.space.Boundary;
import etomica.space.Vector;
import etomica.units.dimensions.*;

/**
 * Calculates orientation-dependent density around the lattice sites of a crystal with radial harmonically mapped averaging
 */
public class MeterDensityAnisotropicRadialHMA implements IDataSource, DataSourceIndependent, AtomLeafAgentManager.AgentSource<Vector> {

    protected final PotentialMaster potentialMaster;
    protected final IteratorDirective id;
    protected DataSourceUniform xDataSourcer;
    protected DataSourceUniform xDataSourcetheta;
    protected DataSourceUniform xDataSourcephi;

    protected final PotentialCalculationForceSum pc;
    protected final AtomLeafAgentManager<Vector> agentManager;
    protected final Box box;
    protected DataFunction data;
    protected IDataInfo dataInfo;
    protected double Rmax;
    protected Vector rivector;

    protected final DataTag tag;
    protected double temperature;
    protected CoordinateDefinition latticesite;
    protected double[] arraymsd;
    protected int rnumberofbins;
    protected int thetaphinumberofbins;

    public MeterDensityAnisotropicRadialHMA(double[] arraymsd, int rnumberofbins, int thetaphinumberofbins, Box box, PotentialMaster potentialMaster, double temperature, CoordinateDefinition latticesite) {
        this.box = box;
        this.temperature = temperature;
        this.arraymsd = arraymsd;
        this.rnumberofbins = rnumberofbins;
        this.thetaphinumberofbins = thetaphinumberofbins;

        this.Rmax = Math.sqrt(arraymsd[0])*4;
        this.latticesite = latticesite;
        this.rivector = box.getSpace().makeVector();
        xDataSourcer = new DataSourceUniform("r", Length.DIMENSION);
        tag = new DataTag();
        xDataSourcer.setTypeMax(LimitType.HALF_STEP);
        xDataSourcer.setTypeMin(LimitType.HALF_STEP);
        xDataSourcer.setXMin(0);
        xDataSourcer.setXMax(Rmax);

        this.potentialMaster = potentialMaster;
        id = new IteratorDirective();
        pc = new PotentialCalculationForceSum();
        agentManager = new AtomLeafAgentManager<Vector>(this, box);
        pc.setAgentManager(agentManager);

        xDataSourcetheta = new DataSourceUniform("theta", Length.DIMENSION);
        xDataSourcetheta.setTypeMax(LimitType.HALF_STEP);
        xDataSourcetheta.setTypeMin(LimitType.HALF_STEP);

        xDataSourcephi = new DataSourceUniform("phi", Length.DIMENSION);
        xDataSourcephi.setTypeMax(LimitType.HALF_STEP);
        xDataSourcephi.setTypeMin(LimitType.HALF_STEP);

        xDataSourcetheta.setXMin(0);
        xDataSourcetheta.setXMax(Math.PI);
        xDataSourcephi.setXMin(0);
        xDataSourcephi.setXMax(2*Math.PI);
        xDataSourcer.setNValues(rnumberofbins);
        xDataSourcetheta.setNValues(thetaphinumberofbins);
        xDataSourcephi.setNValues(thetaphinumberofbins);
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    /**
     * Returns the profile for the current configuration.
     */
    public IData getData() {
        pc.reset();
        potentialMaster.calculate(box, id, pc);
        data.E(0);
        double[] y = data.getData();

        IAtomList atoms = box.getLeafList();
        double dz = Rmax / xDataSourcer.getNValues();
        double[][][] ytemporary = new double[rnumberofbins][thetaphinumberofbins][thetaphinumberofbins];

        for (IAtom atom : atoms) {
                 rivector.Ev1Mv2(atom.getPosition(), latticesite.getLatticePosition(atom));
                 box.getBoundary().nearestImage(rivector);
                 double ri = Math.sqrt(rivector.squared());
                 double fr = agentManager.getAgent(atom).dot(rivector) / ri;
                 double thetai = Math.acos(rivector.getX(2) / ri);
                 double phii = Math.atan2(rivector.getX(1), rivector.getX(0));
                 if (phii < 0) { phii = phii + 2 * Math.PI; }
                 //
                  for (int i = 0; i < rnumberofbins; i++) {
                     double r = (i + 0.5) * dz;

                            int j = (int) (thetai*thetaphinumberofbins/Math.PI);
                     double thetabegin = j * Math.PI / thetaphinumberofbins;
                     double thetaend = (j + 1) * Math.PI / thetaphinumberofbins;

                            int k = (int) (phii*thetaphinumberofbins/(2*Math.PI));
                     double phibegin = k * 2 * Math.PI / thetaphinumberofbins;
                     double phiend = (k + 1) * 2 * Math.PI / thetaphinumberofbins;

                     int mm=j*thetaphinumberofbins+k;

                     double pri = (4.14593*Math.exp(-3 * ri * ri / (2 * arraymsd[mm]))) / ((Math.pow(arraymsd[mm], 1.5))*((phiend-phibegin)*(Math.cos(thetabegin)-Math.cos(thetaend))));
                             double pr = (4.14593*Math.exp(-3 * r * r / (2 * arraymsd[mm]))) / ((Math.pow(arraymsd[mm], 1.5))*((phiend-phibegin)*(Math.cos(thetabegin)-Math.cos(thetaend))));
                             double erf = 1 - SpecialFunctions.erfc(ri * Math.sqrt(3 / (2 * arraymsd[mm])));
                             double cri = (( (ri*1.38198*Math.exp(-3 * ri * ri / (2 * arraymsd[mm])))/(Math.pow(arraymsd[mm], 0.5)) )-(erf))/((phibegin -phiend)*(Math.cos(thetabegin)-Math.cos(thetaend)));
                             double q =1;
                             double dpribydri =-((12.4378*ri*(Math.exp(-3 * ri * ri / (2 * arraymsd[mm]))))/((Math.pow(arraymsd[mm], 2.5))*(phiend-phibegin)*(Math.cos(thetabegin)-Math.cos(thetaend))));


                             double heavisideri;
                             if (ri >= r) { heavisideri = 1; } else { heavisideri = 0; }
                             double diffheavisidethetaend;

                            //UNIFORM MAPPING
                        //     pr=1;pri=1;cri=ri*ri*ri/3;q=500/1.29;dpribydri=0;

                          double ridot = (pr * (((heavisideri ) / ((-phibegin + phiend) * (Math.cos(thetabegin) - Math.cos(thetaend)))) - cri / q)) / (temperature * ri * ri * pri);

                 //           y[n] = y[n] + (pr / q) - (fr - (temperature * dpribydri / pri)) * ridot;
                             ytemporary[i][j][k] = ytemporary[i][j][k] +  (pr / q) - (fr - (temperature * dpribydri / pri)) * ridot;
                 //              System.out.println("n "+n+" "+y[n]+" "+i+" "+j+" "+k );

                 }
             }

        int n=0;
        for (int i = 0; i < rnumberofbins; i++)   {
            for (int j = 0; j < thetaphinumberofbins; j++) {
                for (int k = 0; k < thetaphinumberofbins; k++) {
                    y[n] = ytemporary[i][j][k] ;
                    n=n+1;
                 }
            }
        }

         return data;
    }

    public DataDoubleArray getIndependentData(int i) {
        if(i==0){return (DataDoubleArray) xDataSourcer.getData();}
        if(i==1){return (DataDoubleArray) xDataSourcetheta.getData();}
        if(i==2){return (DataDoubleArray) xDataSourcephi.getData();}
        throw new RuntimeException();
    }
    public DataInfoDoubleArray getIndependentDataInfo(int i) {

        if(i==0){return (DataInfoDoubleArray) xDataSourcer.getDataInfo();}
        if(i==1){return (DataInfoDoubleArray) xDataSourcetheta.getDataInfo();}
        if(i==2){return (DataInfoDoubleArray) xDataSourcephi.getDataInfo();}
        throw new RuntimeException();
    }
    public int getIndependentArrayDimension() {
        return 3;
    }

    public DataTag getIndependentTag(){throw new RuntimeException();}
    public DataTag getIndependentTag(int i) {

        if(i==0){return xDataSourcer.getTag();}
        if(i==1){return xDataSourcetheta.getTag();}
        if(i==2){return xDataSourcephi.getTag();}
        throw new RuntimeException();
    }
    /**
     * @return Returns the box.
     */
    public Box getBox() {
        return box;
    }

    public void reset() {
        if (box == null) return;

        Boundary boundary = box.getBoundary();

        xDataSourcer.setXMin(0);
        xDataSourcer.setXMax(Rmax);
        xDataSourcetheta.setXMin(0);
        xDataSourcetheta.setXMax(Math.PI);
        xDataSourcephi.setXMin(0);
        xDataSourcephi.setXMax(2*Math.PI);

        data = new DataFunction(new int[]{xDataSourcer.getNValues(),xDataSourcetheta.getNValues(),xDataSourcephi.getNValues()});
        dataInfo = new DataInfoFunction("Mapped Average Profile", new CompoundDimension(new Dimension[]{Quantity.DIMENSION, Volume.DIMENSION}, new double[]{1, -1}), this);
        dataInfo.addTag(tag);
    }

    public DataSourceUniform getXDataSource(int i) {
        if(i==0)return xDataSourcer;
        if(i==1)return xDataSourcetheta;
        if(i==2)return xDataSourcephi;
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
