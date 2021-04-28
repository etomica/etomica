/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.molecule.DipoleSourceMolecular;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.TorqueSourceMolecular;
import etomica.potential.PotentialCallbackPhiSumFasterer;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

/**
 * meter for AEE use mapping average
 *
 * @author Weisong
 */
public class MeterDipoleSumSquaredMappedAverageFasterer implements IDataSource {

    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final PotentialCallbackPhiSumFasterer secondDerivativeSum;
    protected final Space space;
    private final Box box;
    private final Vector vectorSum;
    private double dipoleMagnitude;
    private double temperature;
    protected final PotentialCompute potentialMaster;
    protected Vector dr;
    protected Vector work;
    protected TorqueSourceMolecular torqueSourceMolecular;
    protected DipoleSourceMolecular dipoleSource;

    public MeterDipoleSumSquaredMappedAverageFasterer(Box box, double dipoleMagnitude, double temperature,
                                                      PotentialCompute potentialMaster,
                                                      DipoleSourceMolecular dipoleSource) {
        this.dipoleSource = dipoleSource;
        torqueSourceMolecular = new TorqueSourceMolecular(box, potentialMaster);
        data = new DataDoubleArray(2);
        dataInfo = new DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{2});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.box = box;
        this.dipoleMagnitude = dipoleMagnitude;
        this.temperature = temperature;
        this.space = box.getSpace();
        this.potentialMaster = potentialMaster;
        vectorSum = space.makeVector();
//        r = space.makeVector();
        vectorSum.setX(2, 1);
        secondDerivativeSum = new PotentialCallbackPhiSumFasterer(box, dipoleSource);
        dr = space.makeVector();
        work = space.makeVector();

    }

    public IData getData() {        
        double[] x = data.getData();
        double bt = 1/(temperature);
        
        double mu = dipoleMagnitude;
        double mu2 = mu*mu;
        double bt2 = bt*bt;
        double bt3 = bt*bt*bt;

        IMoleculeList moleculeList = box.getMoleculeList();

        int nM = moleculeList.size();
        
        secondDerivativeSum.zeroSum();
        double u = potentialMaster.computeAll(true, secondDerivativeSum);
        if (u == Double.POSITIVE_INFINITY) {
            throw new RuntimeException("u infinity during data collection");
        }

        vectorSum.E(0);
        for (IMolecule m : moleculeList) {
             dr.E(dipoleSource.getDipole(m));
             dr.normalize();

            dr.XE(torqueSourceMolecular.getTorque(m));
            vectorSum.PE(dr);
        }//i loop

        x[0] = -nM*bt2*mu2 - 0.25*bt2*bt2*mu2*vectorSum.squared()+ 0.25*bt3*mu2*secondDerivativeSum.getSum(potentialMaster.getForces());//TODO
        return data;
    }
    
    public DataTag getTag() {
        return tag;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public Box getBox() {
        return box;
    }
}
