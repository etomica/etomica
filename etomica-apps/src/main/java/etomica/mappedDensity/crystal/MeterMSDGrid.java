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
import etomica.space.Vector;
import etomica.units.dimensions.*;

//Meter for calculating 3D MSD

public class MeterMSDGrid implements IDataSource, DataSourceIndependent {


    protected final Box box;
    protected DataSourceUniform xDataSourcetheta;
    protected DataSourceUniform xDataSourcephi;
    protected DataFunction data;
    protected IDataInfo dataInfo;
    protected Vector rivector;
    protected int thetaphinumberofbins;
    protected double ytemporary[][];
    protected int num[][];

    /**
     * Vector describing the orientation of the profile.
     * For example, (1,0) is along the x-axis.
     */
    /**
     * Meter that defines the property being profiled.
     */
    protected final DataTag tag;
    protected CoordinateDefinition latticesite;

    /**
     * Default constructor sets profile along the y-axis, with 100 histogram points.
     */
    public MeterMSDGrid(int thetaphinumberofbins, Box box, CoordinateDefinition latticesite) {
        this.box = box;
        this.thetaphinumberofbins = thetaphinumberofbins;
        this.latticesite = latticesite;
        this.rivector = box.getSpace().makeVector();
        tag = new DataTag();
        ytemporary=new double[thetaphinumberofbins][thetaphinumberofbins];
        num=new int[thetaphinumberofbins][thetaphinumberofbins];

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

        IAtomList atoms = box.getLeafList();
  //      FileWriter fw = null;
  //  try {
  //      fw = new FileWriter("foo.dat",true);
  //  }
  //  catch (IOException e) {throw new RuntimeException(e);}
         for (IAtom atom : atoms) {
            rivector.Ev1Mv2(atom.getPosition(), latticesite.getLatticePosition(atom));
            box.getBoundary().nearestImage(rivector);
             double r2i = rivector.squared();
      //       try {fw.write(""+r2i+"\n");}catch (IOException ex){throw new RuntimeException(ex);}
            double thetai= Math.acos(rivector.getX(2)/Math.sqrt(r2i));
             double phii=Math.atan2(rivector.getX(1),rivector.getX(0));
            if(phii<0){phii=phii+2*Math.PI;}
            int j = (int) (thetai*thetaphinumberofbins/Math.PI);
            int k = (int) (phii*thetaphinumberofbins/(2*Math.PI));
            ytemporary[j][k]=ytemporary[j][k]+r2i;
            num[j][k]++;
         }
  //      try {
  //          fw.close();
  //      }
  //      catch (IOException e) {throw new RuntimeException(e);}

        int n=0;

         for (int j = 0; j < thetaphinumberofbins; j++) {
             for (int k = 0; k < thetaphinumberofbins; k++) {
                 y[n] =ytemporary[j][k]/num[j][k];
                 n=n+1;
              }
         }

         return data;
    }

    public DataDoubleArray getIndependentData(int i) {
        if(i==0){return (DataDoubleArray) xDataSourcetheta.getData();}
        if(i==1){return (DataDoubleArray) xDataSourcephi.getData();}
         throw new RuntimeException();
    }

    public DataInfoDoubleArray getIndependentDataInfo(int i) {

         if(i==0){return (DataInfoDoubleArray) xDataSourcetheta.getDataInfo();}
        if(i==1){return (DataInfoDoubleArray) xDataSourcephi.getDataInfo();}
        throw new RuntimeException();
    }

    public int getIndependentArrayDimension() {
        return 2;
    }

    public DataTag getIndependentTag(int i) {

         if(i==0){return xDataSourcetheta.getTag();}
        if(i==1){return xDataSourcephi.getTag();}
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

         ytemporary=new double[thetaphinumberofbins][thetaphinumberofbins];
        num=new int[thetaphinumberofbins][thetaphinumberofbins];

        xDataSourcetheta.setXMin(0);
        xDataSourcetheta.setXMax(Math.PI);
        xDataSourcephi.setXMin(0);
        xDataSourcephi.setXMax(2*Math.PI);

        data = new DataFunction(new int[]{xDataSourcetheta.getNValues(),xDataSourcephi.getNValues()});
        dataInfo = new DataInfoFunction("Mapped Average Profile", new CompoundDimension(new Dimension[]{Quantity.DIMENSION, Volume.DIMENSION}, new double[]{1, -1}), this);
        dataInfo.addTag(tag);
    }

    public DataSourceUniform getXDataSource(int i) {
         if(i==0)return xDataSourcetheta;
        if(i==1)return xDataSourcephi;
        throw new RuntimeException();

    }

    public DataTag getIndependentTag(){throw new RuntimeException();}

}
