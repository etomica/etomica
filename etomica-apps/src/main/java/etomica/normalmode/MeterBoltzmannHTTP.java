/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.data.histogram.HistogramCollapsing;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.BondingInfo;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.units.dimensions.Null;

/**
 * Same as TargetTP class but is used by MethodBoltzmannHTTP
 * 
 * It is to illustrate the idea of HTTP method.
 * 
 * @author Tai Boon Tan
 *
 */
public class MeterBoltzmannHTTP implements IDataSource {

    protected final PotentialMaster potentialMaster, potentialMasterPretend;
    protected double latticeEnergy;
    protected double temperature;
    protected double otherTemperature;
    protected DataDoubleArray data;
    protected DataInfoDoubleArray dataInfo;

    protected DataTag tag;

    protected final Box pretendBox;
    protected final CoordinateDefinition coordinateDefinition;
    protected final ISpecies species;
    protected P1ConstraintNbr p1;
    protected HistogramCollapsing[] histogram;

    public MeterBoltzmannHTTP(PotentialMaster potentialMaster, ISpecies species, Simulation sim, CoordinateDefinition coordinateDefinition, double rc) {
        this.potentialMaster = potentialMaster;
        this.coordinateDefinition = coordinateDefinition;
        this.species = species;
        pretendBox = sim.makeBox(coordinateDefinition.getBox().getBoundary());
        potentialMasterPretend = new PotentialMasterList(sim.getSpeciesManager(),
                pretendBox, 7, rc, BondingInfo.noBonding(), potentialMaster.getPairPotentials());

        data = new DataDoubleArray(5);
        dataInfo = new DataInfoDoubleArray("Reduced Energy", Null.DIMENSION, new int[]{5});

        tag = new DataTag();

        // insert atoms into the box at their lattice sites.
        // we do this because want to find neighbors now (and then never again)
        Box realBox = coordinateDefinition.getBox();
        pretendBox.setNMolecules(species, realBox.getNMolecules(species));
        IAtomList atoms = realBox.getLeafList();
        IAtomList pretendAtoms = pretendBox.getLeafList();
        for (int j = 0; j<atoms.size(); j++) {
            IAtom jRealAtom = atoms.get(j);
            Vector pos = pretendAtoms.get(j).getPosition();
            pos.E(coordinateDefinition.getLatticePosition(jRealAtom));
        }
        potentialMasterPretend.init();
    }
    
    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    public IData getData() {
        Box realBox = coordinateDefinition.getBox();
        double u = potentialMaster.computeAll(false);

        IAtomList atoms = realBox.getLeafList();
        IAtomList pretendAtoms = pretendBox.getLeafList();
        double a0 = (u-latticeEnergy)/temperature;
        double a1NotScaled = (u-latticeEnergy)/otherTemperature;
        
        double[] x = data.getData();

        double fac = Math.sqrt(otherTemperature/temperature);
        
        for (int j = 0; j<atoms.size(); j++) {
            IAtom jRealAtom = atoms.get(j);
            Vector pos = pretendAtoms.get(j).getPosition();
            pos.Ea1Tv1(1-fac, coordinateDefinition.getLatticePosition(jRealAtom));
            pos.PEa1Tv1(+fac, jRealAtom.getPosition());
        }

        double otherU = 0;
        if (p1 != null) {
            // we need to include the constraint energy here even though we
            // didn't include it in the u above (for realBox).  RealBox will
            // always have constraint energy = 0 (moves that violate the
            // constraint are rejected).  But we could have scaled the
            // atoms into a constraint violation, so we should check now.
            otherU = constraintEnergy(pretendBox);
        }
        if (otherU < Double.POSITIVE_INFINITY) {
            otherU += potentialMasterPretend.computeAll(false);
        }
        double a1 = (otherU-latticeEnergy)/otherTemperature;
        
        x[0] = a0;
        x[1] = a1NotScaled;
        x[2] = a1;
        
        x[3] = a1NotScaled-a0;
        x[4] = a1-a0;	
        
        return data;
    }

    /**
     * Returns true if all atoms in the given box satisfy p1's constraint
     */
    protected double constraintEnergy(Box box) {
        p1.setBox(box);
        IAtomList atomList = box.getLeafList();
        for (int i = 0; i<atomList.size(); i++) {
            if (p1.energyi(atomList.get(i)) == Double.POSITIVE_INFINITY) {
                return Double.POSITIVE_INFINITY;
            }
        }
        // ugly!!!
        p1.setBox(coordinateDefinition.getBox());
        return 0;
    }

   
    public double getLatticeEnergy() {
        return latticeEnergy;
    }

    public void setLatticeEnergy(double latticeEnergy) {
        this.latticeEnergy = latticeEnergy;
    }

    public double getTemperature() {
        return temperature;
    }

    public void setTemperature(double temperature) {
        this.temperature = temperature;
    }

    public double getOtherTemperature() {
        return otherTemperature;
    }

    public void setOtherTemperature(double otherTemperature) {
        this.otherTemperature = otherTemperature;
    }
    
    public CoordinateDefinition getCoordinateDefinition() {
        return coordinateDefinition;
    }

    public void setConstraint(P1ConstraintNbr p1) {
        this.p1 = p1;
        p1.initBox(pretendBox);
    }
}
