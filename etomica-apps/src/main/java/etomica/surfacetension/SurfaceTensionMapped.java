/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.surfacetension;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataProcessor;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.meter.MeterNMolecules;
import etomica.data.meter.MeterProfileByVolume;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataDoubleArray;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.units.dimensions.Null;


/**
 * DataProcessor that adds in mapped average contributions to the surface
 * tension.
 * 
 * @author Andrew Schultz
 */
public class SurfaceTensionMapped extends DataProcessor {

    protected final DataDouble data = new DataDouble();
    protected final FitTanh fit;
    protected final MeterProfileByVolume densityProfileMeter;
    protected final PotentialMaster potentialMaster;
    protected final PotentialCalculationForceSum pcForce;
    protected final AtomLeafAgentManager<Vector> forceAgentManager;
    protected final Space space;
    protected final Box box;
    protected final IteratorDirective allAtoms;

    public SurfaceTensionMapped(Space space, Box box, ISpecies species, PotentialMaster potentialMaster) {
        this.space = space;
        this.box = box;
        densityProfileMeter = new MeterProfileByVolume(space);
        densityProfileMeter.setBox(box);
        densityProfileMeter.getXDataSource().setNValues(400);
        MeterNMolecules meterNMolecules = new MeterNMolecules();
        meterNMolecules.setSpecies(species);
        densityProfileMeter.setDataSource(meterNMolecules);

        fit = new FitTanh();
        
        this.potentialMaster = potentialMaster;
        pcForce = new PotentialCalculationForceSum();
        forceAgentManager = new AtomLeafAgentManager<>(a -> space.makeVector(), box);
        pcForce.setAgentManager(forceAgentManager);
        allAtoms = new IteratorDirective();
    }

    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        dataInfo = new DataInfoDouble("surface tension", Null.DIMENSION);
        return dataInfo;
    }

    protected IData processData(IData inputData) {
        double st = inputData.getValue(0);
        
        double[] x = ((DataDoubleArray)densityProfileMeter.getXDataSource().getData()).getData();
        double[] y = ((DataDoubleArray)densityProfileMeter.getData()).getData();
        double[] param = fit.doFit(x, y);
        double c = 2/param[2];
//        System.out.println("fit: "+Arrays.toString(param));
        IAtomList atoms = box.getLeafList();
        double rho = atoms.size()/box.getBoundary().volume();
        double L = box.getBoundary().getBoxSize().getX(0);
        double tL1 = Math.tanh(c * (L/2 - param[3]));
        double tL2 = Math.tanh(c * (L/2 + param[3]));
        double rhoL2 = param[0] + 0.5 * (param[1] - param[0]) * (tL2 - tL1);
        double dzsdL = (rho-rhoL2)/(param[1]-param[0])/(tL2 + tL1);

        pcForce.reset();;
        potentialMaster.calculate(box, allAtoms, pcForce);
        double mapSum = 0;



        double cL1 = Math.cosh(c*(L/2-param[3]));
        double cL2 = Math.cosh(c*(L/2+param[3]));
        double jFac = (param[1]-param[0])/(2*L) * ((tL2-tL1)*c*L + 2*Math.log(cL1/cL2))/(tL1+tL2);
        for (int i = 0; i<atoms.size(); i++) {
            IAtom atom = atoms.get(i);
            Vector f = forceAgentManager.getAgent(atom);
            double px = atom.getPosition().getX(0);
            double t1 = Math.tanh(c * (px + param[3]));
            double t2 = Math.tanh(c * (px - param[3]));
            double irho = param[0] + 0.5 * (param[1] - param[0]) * (t1 - t2);
            double cosh1 = Math.cosh(c*(px-param[3]));
            double cosh2 = Math.cosh(c*(px+param[3]));
            double foo = param[0]*px/L - (param[1]-param[0])/(2*c*L) * Math.log(cosh1/cosh2);
            double bar = dzsdL * 0.5 * (param[1]-param[0]) * (t1 + t2);
            double xL = (foo - bar)/irho; 
            mapSum += (xL - px/L)*f.getX(0);

            double sech1 = 1.0/Math.cosh(c*(px-param[3]));
            double sech2 = 1.0/Math.cosh(c*(px+param[3]));
            mapSum += jFac*(sech1*sech1+sech2*sech2)/(2*rho + (rho-param[1])*(t2 - t1));
        }
        
        data.x = st + mapSum/box.getBoundary().volume();
        return data;
    }

}
