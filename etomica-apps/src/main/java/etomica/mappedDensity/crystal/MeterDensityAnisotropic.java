/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.mappedDensity.crystal;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.DataSourceUniform.LimitType;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.normalmode.CoordinateDefinition;
import etomica.space.Boundary;
import etomica.space.Vector;
import etomica.units.dimensions.*;

/**
 *
 * Histograms orientation-dependent density around the lattice sites of a crystal
 */

public class MeterDensityAnisotropic implements IDataSource, DataSourceIndependent {


    protected final Box box;
    protected DataSourceUniform xDataSourcer;
    protected DataSourceUniform xDataSourcetheta;
    protected DataSourceUniform xDataSourcephi;

    protected DataFunction data;
    protected IDataInfo dataInfo;
    protected double Rmax;
    protected Vector rivector;
    protected int rnumberofbins;
    protected int thetaphinumberofbins;

    protected final DataTag tag;
    protected CoordinateDefinition latticesite;

    public MeterDensityAnisotropic(double[] arraymsd, int rnumberofbins, int thetaphinumberofbins, Box box, CoordinateDefinition latticesite) {
        this.box = box;
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
//sim.coordinateDefinition - latticesites
    /**
     * Accessor method for vector describing the direction along which the profile is measured.
     * Each atom position is dotted along this vector to obtain its profile abscissa value.
     */
    /**
     * Accessor method for vector describing the direction along which the profile is measured.
     * Each atom position is dotted along this vector to obtain its profile abscissa value.
     * The given vector is converted to a unit vector, if not already.
     */
    //       double pminusLby2= c.df(1, -L/2);
    //       double zidotminusLby2=-0.0725;
    //     zidotminusLby2=pz*(cz/q-1/2)/(pminusLby2*temperature);
    //      if (zi >= z) {zidotminusLby2=-pz*(1-(czplusLby2/q))/(pminusLby2*temperature);} else {zidotminusLby2=pz*((czminusLby2/q))/(pminusLby2*temperature);}
    //      double heavisidei;
    //      if (zi >= z) {heavisidei=1;} else {heavisidei=0;}
    //      return ((pminusLby2*zidotminusLby2/pzi)+(pz*((heavisidei)-(czi/q))/(temperature*pzi)));
    //    return ((((heavisidei)-(zi/L))/(temperature)));


    /**
     * Returns the profile for the current configuration.
     */
    public IData getData() {
        data.E(0);
        double[] y = data.getData();
        double [][][] ytemporary=new double[rnumberofbins][thetaphinumberofbins][thetaphinumberofbins];

    //    for (int i = 0; i < rnumberofbins; i++) {
    //        for (int j = 0; j < thetaphinumberofbins; j++) {
    //            for (int k = 0; k < thetaphinumberofbins; k++) {
    //                ytemporary[i][j][k]=0;
    //            }
    //        }
    //    }
            IAtomList atoms = box.getLeafList();
        double dz = Rmax / xDataSourcer.getNValues();

         for (IAtom atom : atoms) {
            rivector.Ev1Mv2(atom.getPosition(), latticesite.getLatticePosition(atom));
            box.getBoundary().nearestImage(rivector);
            double ri = Math.sqrt(rivector.squared());
            double thetai= Math.acos(rivector.getX(2)/ri);
             double phii=Math.atan2(rivector.getX(1),rivector.getX(0));
            if(phii<0){phii=phii+2*Math.PI;}
            int i = (int) (ri/dz);
            int j = (int) (thetai*thetaphinumberofbins/Math.PI);
            int k = (int) (phii*thetaphinumberofbins/(2*Math.PI));
            ytemporary[i][j][k]=ytemporary[i][j][k]+1;
        }

        int n=0;
        int numAtoms = box.getLeafList().size();
     for (int i = 0; i < rnumberofbins; i++)   {
            double Rbinstart=i*dz;
         double Rbinend=(i+1)*dz;

         for (int j = 0; j < thetaphinumberofbins; j++) {
             double thetabegin = j * Math.PI / thetaphinumberofbins;
             double thetaend = (j + 1) * Math.PI / thetaphinumberofbins;

             for (int k = 0; k < thetaphinumberofbins; k++) {
                 double phibegin = k * 2 * Math.PI / thetaphinumberofbins;
                 double phiend = (k + 1) * 2 * Math.PI / thetaphinumberofbins;

                 y[n] = ytemporary[i][j][k] / numAtoms / (((-phibegin + phiend) * (Math.cos(thetabegin) - Math.cos(thetaend))) * (Rbinend * Rbinend * Rbinend - Rbinstart * Rbinstart * Rbinstart) / 3);
                 n=n+1;

                 //   System.out.println("qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq "+q);
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
 //       xDataSource.setYMin(0);

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

    public DataTag getIndependentTag(){throw new RuntimeException();}

}
