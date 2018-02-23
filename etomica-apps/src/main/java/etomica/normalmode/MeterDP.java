/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.math.function.Function;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.units.dimensions.Null;

import java.io.FileWriter;

/**
 * Meter that measures the overlap averages for perturbing from a solid at one
 * density into other densities.  The atoms are scaled back toward their
 * lattice sites by a factor that depends on the density change and an estimate
 * deltaA/deltaV.  Different scalings are applied for different estimates of
 * deltaA/deltaV, with alpha=1 used (for the overlap average) for each scaling.
 */
public class MeterDP implements IDataSource {

    protected final MeterPotentialEnergy meterPotential;
    protected final PotentialMaster potentialMaster;
    protected double temperature;
    protected double[] otherDensities;
    protected DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected final DataTag tag;
    protected final Box pretendBox;
    protected CoordinateDefinition coordinateDefinition;
    protected final ISpecies species;
    protected double[][] p;
    protected double[] pCenter;
    protected double pSpan;
    protected boolean linPSpan;
    protected int numP = 1;
    protected P1ConstraintNbr p1;
    protected final Vector work;
    protected Function uLatFunction = uLat0;
    protected FileWriter fw;
    
    public MeterDP(PotentialMaster potentialMaster, ISpecies species, Space space, Simulation sim) {
        this.potentialMaster = potentialMaster;
        meterPotential = new MeterPotentialEnergy(potentialMaster);
        this.species = species;
        pretendBox = sim.makeBox();
        if (potentialMaster instanceof PotentialMasterList) {
            pretendBox.getBoundary().getEventManager().removeListener(((PotentialMasterList) potentialMaster).getNbrCellManager(pretendBox));
        }

        work = space.makeVector();

        tag = new DataTag();
//        try {
//            fw = new FileWriter("dp.dat");
//        }
//        catch (IOException e) {throw new RuntimeException(e);}
    }
    
    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    /**
     * Sets a function that returns the lattice energy for a given density.
     * If not set, the default lattice energy is taken to be 0 for all densities.
     */
    public void setULatFunction(Function newULatFunction) {
        uLatFunction = newULatFunction;
    }
    
    public Function getULatFunction() {
        return uLatFunction;
    }
    
    public IData getData() {
        Box realBox = coordinateDefinition.getBox();
        meterPotential.setBox(realBox);
        double u = 0; //meterPotential.getDataAsScalar();
        
        int nMolecules = realBox.getMoleculeList().getMoleculeCount();
        double rho = nMolecules / realBox.getBoundary().volume();
        int D = realBox.getBoundary().getBoxSize().getD();
        
        meterPotential.setBox(pretendBox);

        IAtomList atoms = realBox.getLeafList();
        IAtomList pretendAtoms = pretendBox.getLeafList();
        int numAtoms = atoms.size();
        double uLatRho0 = numAtoms*uLatFunction.f(rho);
        double a0 = (u-uLatRho0)/temperature;

        double[] x = data.getData();
        double v0 = realBox.getBoundary().volume();
        for (int i=0; i<otherDensities.length; i++) {
            double rScale = Math.pow(rho / otherDensities[i], 1.0/D);
            work.Ea1Tv1(rScale, realBox.getBoundary().getBoxSize());
            pretendBox.getBoundary().setBoxSize(work);
            double vi = pretendBox.getBoundary().volume();
            double uLatRhoi = numAtoms*uLatFunction.f(otherDensities[i]);

            for (int k=0; k<numP; k++) {
                double fac = Math.exp((p[i][k]*(vi-v0) + (uLatRhoi-uLatRho0))/(numAtoms*temperature*D));
                if (fac > 2) continue;
//                System.out.println(i+" "+k+" "+p[i][k]+" "+rScale+" "+fac+" "+((p[i][k]*(vi-v0) + (uLatRhoi-uLatRho0))/(numAtoms*temperature*D)));
                for (int j=0; j<numAtoms; j++) {
                    IAtom jRealAtom = atoms.get(j);
                    Vector pos = pretendAtoms.get(j).getPosition();
//                    if (j==25) System.out.println("hi "+coordinateDefinition.getLatticePosition(jRealAtom)+" "+jRealAtom.getPosition());
                    pos.Ea1Tv1(rScale-fac, coordinateDefinition.getLatticePosition(jRealAtom));
                    pos.PEa1Tv1(+fac, jRealAtom.getPosition());
//                    if (j==25) System.out.println("newP "+pos);
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
                    otherU += meterPotential.getDataAsScalar();
                }
                double ai = (otherU-uLatRhoi)/temperature;
//                System.out.println(p[i][k]+" "+u+" "+uLatRho0+" "+otherU+" "+uLatRhoi+" "+a0+" "+ai);
                x[i*numP+k] = 1.0/(1+Math.exp(ai-a0));
//                if (i==0 && k==(numP-1)/2) {
//                    try {
//                        fw.write(x[i*numP+k]+"\n");
//                    }
//                    catch (IOException e) {throw new RuntimeException(e);}
//                }
            }
        }
//        System.out.println(x[1]);
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
        return 0;
    }

    public double getTemperature() {
        return temperature;
    }

    public void setTemperature(double temperature) {
        this.temperature = temperature;
    }

    public double[] getOtherDensities() {
        return otherDensities;
    }

    public void setOtherDensities(double[] otherDensities) {
        this.otherDensities = otherDensities;
    }
    
    protected void initP() {
        if (pCenter == null) {
            return;
        }
        p = new double[pCenter.length][numP];
        if (numP == 1) {
            for (int i=0; i<pCenter.length; i++) {
                p[i][0] = pCenter[i];
            }
        }
        else {
            for (int i=0; i<pCenter.length; i++) {
                for (int j=0; j<numP; j++) {
                    if (linPSpan) {
                        p[i][j] = pCenter[i] + 2.0*pSpan*(j-(numP-1)/2)/(numP-1);
                    }
                    else {
                        p[i][j] = pCenter[i] * Math.exp(2.0*pSpan*(j-(numP-1)/2)/(numP-1));
                    }
                }
            }
        }
        data = new DataDoubleArray(numP*pCenter.length);
        dataInfo = new DataInfoDoubleArray("overlap", Null.DIMENSION, new int[]{numP*pCenter.length});
    }

    public double[] getP(int iP) {
        return p[iP];
    }
    
    public double getPCenter(int iTemp) {
        return pCenter[iTemp];
    }
    
    public void setPCenter(double[] newPCenter) {
        pCenter = newPCenter;
    }
    
    public void setExpPSpan(double newPSpan) {
        pSpan = newPSpan;
        linPSpan = false;
        initP();
    }
    
    public void setLinPSpan(double newPSpan) {
        pSpan = newPSpan;
        linPSpan = true;
        initP();
    }
    
    public void setNumP(int newNumP) {
        numP = newNumP;
        initP();
    }

    public CoordinateDefinition getCoordinateDefinition() {
        return coordinateDefinition;
    }

    public void setCoordinateDefinition(CoordinateDefinition newCoordinateDefinition) {
        this.coordinateDefinition = newCoordinateDefinition;

        // insert atoms into the box at their lattice sites.
        // we do this because want to find neighbors now (and then never again)
        Box realBox = coordinateDefinition.getBox();
        pretendBox.getBoundary().setBoxSize(realBox.getBoundary().getBoxSize());
        pretendBox.setNMolecules(species, realBox.getNMolecules(species));
        IAtomList atoms = realBox.getLeafList();
        IAtomList pretendAtoms = pretendBox.getLeafList();
        for (int j = 0; j<atoms.size(); j++) {
            IAtom jRealAtom = atoms.get(j);
            Vector pos = pretendAtoms.get(j).getPosition();
            pos.E(coordinateDefinition.getLatticePosition(jRealAtom));
        }

        if (potentialMaster instanceof PotentialMasterList) {
            // find neighbors now.
            ((PotentialMasterList)potentialMaster).getNeighborManager(pretendBox).reset();
        }
        
        if (p1 != null) {
            p1.initBox(pretendBox);
        }
    }

    public void setConstraint(P1ConstraintNbr p1) {
        this.p1 = p1;
        if (coordinateDefinition != null) {
            p1.initBox(pretendBox);
        }
    }

    /**
     * Nominal function for lattice energy
     */
    public final static Function uLat0 = new Function.Constant(0);
}
