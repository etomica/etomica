/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.IRandom;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.atom.AtomHydrogen;
import etomica.atom.AtomTypeOrientedSphere;
import etomica.atom.IAtomOriented;
import etomica.atom.IAtomTypeOriented;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.chem.elements.Hydrogen;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.simulation.Simulation;
import etomica.space.ISpace;
import etomica.space.IVectorRandom;
import etomica.space3d.IOrientation3D;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresHetero;
import etomica.units.Kelvin;
import etomica.util.Constants;
import etomica.util.DoubleRange;
import etomica.util.HistogramNotSoSimple;
import etomica.util.HistogramSimple;

/**
 * MCMove that fully regrows the beads of a ring polymer, accepting or
 * rejecting the move based on the sampling weight.  The move can (optionally)
 * regrow the beads such that the beads from multiple molecules are combined
 * to form a larger ring.
 * 
 * @author Andrew Schultz
 */
public class MCMoveClusterRingRegrowOrientation extends MCMoveBox {

    
    public MCMoveClusterRingRegrowOrientation(IRandom random, ISpace _space, int P) {    
        super(null);
        this.space = _space;
        this.P = P;
        this.random = random;                
        leafIterator = new AtomIteratorLeafAtoms();
        dr = space.makeVector();
        dr1 = space.makeVector();
        dr2 = space.makeVector();
        dummy = space.makeVector();
	}  
    public void setBox(IBox p) {
        super.setBox(p);
        int nMolecules = box.getMoleculeList().getMoleculeCount();
        xOld = new double[nMolecules][P];
        xNew = new double[nMolecules][P];
        oldOrientations = new IOrientation3D[nMolecules][];
        for (int i=0; i<nMolecules; i++) {
            int nAtoms = box.getMoleculeList().getMolecule(i).getChildList().getAtomCount();
            oldOrientations[i] = new IOrientation3D[nAtoms+1];
            for (int j=0; j<nAtoms+1; j++) {
                oldOrientations[i][j] = (IOrientation3D) space.makeOrientation();
            }
        }
        leafIterator.setBox(p);
//        System.out.println(1E-24*Constants.AVOGADRO+" "+Mole.UNIT.fromSim(1)+" "+1/Constants.AVOGADRO);
//        System.exit(1);
    }
    
    /**
     * Set PI harmonic spring stiffness for given temperature and atomic mass.
     * Dimer is assumed to be composed of two atoms of the given mass.
     */
    public void setStiffness(double t, double mass) {
        double lambda = Constants.PLANCK_H/(Math.sqrt(2*Math.PI*mass*t));
        stiffness = Math.PI*P/(lambda*lambda);        
//        System.out.println("Stiffness = "+stiffness);
    }
    public double getStiffness() {
        return stiffness;
    }
	public boolean doTrial() {
//	    if (oldFancy == null && printToFile) {            
//            try {                
//                oldFancy = new FileWriter("oldFancy.dat");
//                newFancy = new FileWriter("newFancy.dat");
//            } catch (IOException e) {
//                throw new RuntimeException(e);
//            }
//        }
        weightOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);        
        IVectorRandom e1 = (IVectorRandom) space.makeVector();
        IVectorMutable ex = space.makeVector();
        IVectorMutable ey = space.makeVector();
        IVectorMutable ez = space.makeVector();
        ex.setX(0, 1);
        ey.setX(1, 1);
        ez.setX(2, 1);
        IVectorMutable oldCenter = space.makeVector();
        IVectorMutable newCenter = space.makeVector();                
        IMoleculeList molecules = box.getMoleculeList();
        IOrientation3D [][] newOrientations = new IOrientation3D[molecules.getMoleculeCount()][P+1];                
        double [] oldAlpha = new double [P];
        newAlpha = new double [P];
        double [] theta = new double [P];        
        int fromImage = 0;
        int toImage = 0;
        double uOld = 0;
        double uNew = 0;
        double pGenOld = 1.00;
        double pGenNew = 1.00;
        double pGenRatio = 1.00;
        IVectorMutable pVecOld = space.makeVector();
        IVectorMutable pVecNew = space.makeVector();
        IVectorMutable prevOldCenter = space.makeVector();        
        
        for (int i=0; i<molecules.getMoleculeCount(); i++) {
            double bondLength = 0.00;
            IMolecule molecule = molecules.getMolecule(i);
            IAtomList atoms = molecule.getChildList();
            for (int j=0; j<P; j++) {                                
                int prev = j-1;
                if (prev < 0) prev = P-1;
                AtomHydrogen jAtom = (AtomHydrogen)atoms.getAtom(j);
                AtomHydrogen jPrev = (AtomHydrogen)atoms.getAtom(prev);
                uOld += stiffness*dist(jAtom,jPrev);
            }
            
//            System.out.println(BohrRadius.UNIT.fromSim(bondLength));
//            double phi02old = Math.acos(((AtomHydrogen)atoms.getAtom(0)).getOrientation().getDirection().dot(((AtomHydrogen)atoms.getAtom(2)).getOrientation().getDirection()));
//            if (printToFile) {
//                try {
//                    oldFancy.write(phi02old+" ");
//                }
//                catch (IOException e) {
//                    throw new RuntimeException(e);
//                }
//            }
            oldOrientations[i][0].setDirection(((IAtomOriented) atoms.getAtom(0)).getOrientation().getDirection());
            IVectorRandom rV1 = (IVectorRandom)space.makeVector();            
            rV1.setRandomSphere(random);
//            rV1.E(ex);
            newOrientations[i][0] = (IOrientation3D)((IAtomOriented) atoms.getAtom(0)).getOrientation();
            newOrientations[i][0].setDirection(rV1);
            oldOrientations[i][P].setDirection(oldOrientations[i][0].getDirection());
            newOrientations[i][P] = (IOrientation3D) space.makeOrientation();
            newOrientations[i][P].setDirection(newOrientations[i][0].getDirection());
            pVecOld.E(oldOrientations[i][0].getDirection());
            pVecNew.E(newOrientations[i][0].getDirection());
            for (int dr = 2; dr<= P; dr*=2){
                double kEff = 8*stiffness*dr/P;
                double y0 = 0;
                double sA = 0;
                double kEff1Old = 0;
                double kEff1New = 0;
                double a = 0;
                double y1 = 0;
                for (int nr = 1; nr<dr; nr+=2){
                    int imageIndex = nr*P/dr;
//                    System.out.println("image # = "+imageIndex);                    
                    IAtomOriented jAtom = ((IAtomOriented)atoms.getAtom(imageIndex));                   
                    oldOrientations[i][imageIndex].setDirection(jAtom.getOrientation().getDirection());
                    newOrientations[i][imageIndex] = (IOrientation3D) jAtom.getOrientation();
                    fromImage = (nr-1)*P/dr;
                    toImage = (nr+1)*P/dr;
                    double r = 0;
                    for (int k = fromImage; k<= toImage; k++) {
                        if (k == P) {
                            r += ((AtomHydrogen)atoms.getAtom(0)).getBondLength()/2.0;                            
                        }
                        else {
                            r += ((AtomHydrogen)atoms.getAtom(k)).getBondLength()/2.0;
                        }
                        
                    }
                    r *= dr/(2.0*P + dr); // same as r /= (toImage - fromImage + 1)
                    
                    if (imageIndex == P/2 && doExchange) {
                    	pVecOld.TE(-1);
                    	pVecNew.TE(-1);
                    	oldOrientations[i][P].setDirection(pVecOld);
                    	newOrientations[i][P].setDirection(pVecNew);
                    }                    
                    else {                    
                    	oldCenter.Ev1Pv2(oldOrientations[i][fromImage].getDirection(), oldOrientations[i][toImage].getDirection());                    
                    	oldCenter.normalize();
                    	y0 = oldOrientations[i][fromImage].getDirection().dot(oldOrientations[i][toImage].getDirection());
                    	if (y0 > 1.0) y0 = 1.0;
                    	if (y0 < -1.0) y0 = -1.0;
                    	kEff1Old = kEff*r*r*Math.sqrt((1+y0)/2.0);                    
                    	y0 = newOrientations[i][fromImage].getDirection().dot(newOrientations[i][toImage].getDirection());
                    	if (y0 > 1.0) y0 = 1.0;
                    	if (y0 < -1.0) y0 = -1.0;
                    	kEff1New = kEff*r*r*Math.sqrt((1+y0)/2.0);                    
                    	y0 = oldCenter.dot(oldOrientations[i][imageIndex].getDirection());                    
                    	if (y0 > 1.0) y0 = 1.0;
                    	if (y0 < -1.0) y0 = -1.0;
                    	oldAlpha[imageIndex] = Math.acos(y0);                    
                    	double x = random.nextDouble();
                    	xNew[i][imageIndex] = x;                    
//                  	double y1_new = (-kEff1New + Math.log(Math.exp(2*kEff1New)-x*(Math.exp(2*kEff1New)-1)))/kEff1New;                                        
                    	a = Math.log(1 - x)/kEff1New + Math.log(1 + x*Math.exp(-2*kEff1New)/(1-x))/kEff1New;                    
//                  	double a = Math.log(1 - xNew[i][imageIndex])/kEff1New + Math.log(1 + xNew[i][imageIndex]*Math.exp(-2*kEff1New)/(1-xNew[i][imageIndex]))/kEff1New;
//                  	if (firstMove && imageIndex == 2) y1 = Math.cos(Degree.UNIT.toSim(20));                                        
                    	if (a > 0) {                        
                    		a = 0;
                    	}
                    	if (a < -2.0) {                        
                    		a = -2.0;
                    	}
                    	y1 = 1 + a;
                    	sA = Math.sqrt(-2*a - a*a);
                    	if (Double.isNaN(sA)) throw new RuntimeException(a+" "+(2*a + a*a));
                    	newAlpha[imageIndex] = Math.acos(y1);
                    	if (newAlpha[imageIndex] != newAlpha[imageIndex] || y1 != y1) throw new RuntimeException("x = " +2*kEff1New);
//                  	if (imageIndex == 2) hist1.addValue(newAlpha[imageIndex]);
                    }
                    if (imageIndex == P/2) {
                        newOrientations[i][imageIndex].setDirection(rV1);
                        IVectorMutable rV2 = space.makeVector();
                        if (Math.abs(rV1.getX(0)) > 0.5) {
                            rV2.setX(1,1);
                        }
                        else {
                            rV2.setX(0, 1);
                        }
                        rV2.PEa1Tv1(-rV2.dot(rV1), rV1);
                        rV2.normalize();
                        double dummyAlpha = 2*Math.PI*random.nextDouble();                        
                        rotateBy(Math.cos(dummyAlpha), Math.sin(dummyAlpha), rV1, rV2);
//                        newAlpha[imageIndex] = Degree.UNIT.toSim(0);
//                        y1 = 1.0;                        
                        if (!doExchange) {
                        	rotateBy(y1,sA,rV2,(IVectorMutable)newOrientations[i][imageIndex].getDirection());                        
                        }
                        else {
                        	double angle = 2*Math.PI*random.nextDouble();
                        	rotateBy(Math.cos(angle), Math.sin(angle),rV2,(IVectorMutable)newOrientations[i][imageIndex].getDirection());
                        }
                        theta[imageIndex] = 0;
                                                
                    }
                    else {
                        newCenter.Ev1Pv2(newOrientations[i][fromImage].getDirection(), newOrientations[i][toImage].getDirection());
                        newCenter.normalize();
                        newOrientations[i][imageIndex].setDirection(newCenter);                        
//                        System.out.println(newOrientations[i][imageIndex].getDirection());
                        e1.E(0);
                        if (Math.abs(newCenter.getX(0)) > 0.5) {
                            e1.setX(1,1);
                        }
                        else {
                            e1.setX(0, 1);
                        }
                        e1.PEa1Tv1(-e1.dot(newCenter), newCenter);
                        e1.normalize();
//                        double alpha = 2*Math.PI*random.nextDouble();                        
//                        rotateBy(Math.cos(alpha), Math.sin(alpha),newCenter, e1);
                        
                        rotateBy(y1,sA, e1,(IVectorMutable)newOrientations[i][imageIndex].getDirection());
                        theta[imageIndex] = random.nextDouble()*2*Math.PI;                
                        newOrientations[i][imageIndex].rotateBy(theta[imageIndex], newCenter);
//                        double phi02new = Math.acos(((AtomHydrogen)atoms.getAtom(0)).getOrientation().getDirection().dot(((AtomHydrogen)atoms.getAtom(2)).getOrientation().getDirection()));
//                        if (Math.abs(phi02new - newAlpha[2]) > 1.0E-10 ) {
//                            System.out.println(phi02new+" "+newAlpha[2]+" "+e1+" "+newCenter);
//                            throw new RuntimeException();
//                        }
//                      
                    }
                    if (newOrientations[i][imageIndex].getDirection().isNaN()) throw new RuntimeException("bead "+imageIndex+" orientation is NaN");
                    if (!doExchange || imageIndex != P/2) {
                    	oldCenter.ME(oldOrientations[i][imageIndex].getDirection());
                    	double v = oldCenter.squared();
                    	double v1 = oldCenter.dot(oldOrientations[i][imageIndex].getDirection());
                    	double s1 = v - v1*v1;
                    	double y0m1 = -s1/(1+y0);
//                  	if (xOld[i][imageIndex] == 0) xOld[i][imageIndex] = xNew[i][imageIndex];
//                  	pGenOld *= Math.exp(kEff1Old*y0)*kEff1Old/Math.sinh(kEff1Old);                    
//                  	pGenNew *= Math.exp(kEff1New*y1)*kEff1New/Math.sinh(kEff1New);
//                  	pGenOld *= 2*Math.exp(kEff1Old*y0m1)*kEff1Old/(1 - Math.exp(-2*kEff1Old));                    
//                  	pGenNew *= 2*Math.exp(kEff1New*a)*kEff1New/(1 - Math.exp(-2*kEff1New));
//                  	System.out.println(kEff1Old+" "+kEff1New+" "+y0m1+" "+(y0-1)+" "+a+" "+(y1_new-1));
//                  	pGenRatio *= Math.exp(kEff1New*a - kEff1Old*y0m1)*kEff1New*(1-Math.exp(-2*kEff1Old))/(kEff1Old*(1-Math.exp(-2*kEff1New)));
//                  	double aNew = kEff1New*(-xNew[i][imageIndex] + 1/(1 - Math.exp(-2*kEff1New)))/(kEff1Old*(-xOld[i][imageIndex] + 1/(1 - Math.exp(-2*kEff1Old))));
                    	double aOld =  Math.exp(kEff1New*a - kEff1Old*y0m1)*kEff1New*(1-Math.exp(-2*kEff1Old))/(kEff1Old*(1-Math.exp(-2*kEff1New)));                    
//                  	System.out.println(aOld+" "+aNew);
                    	pGenRatio *= aOld;                    
                    	if (Double.isNaN(pGenRatio)) throw new RuntimeException();
                    }
                }
            }            
            for (int j=0; j<P; j++) {
                int prev = j-1;
                if (prev < 0) prev = P-1;
                AtomHydrogen jAtom = (AtomHydrogen)atoms.getAtom(j);
                AtomHydrogen jPrev = (AtomHydrogen)atoms.getAtom(prev);
                uNew += stiffness*dist(jAtom,jPrev);
            }
            
//            double phi02new = Math.acos(((AtomHydrogen)atoms.getAtom(0)).getOrientation().getDirection().dot(((AtomHydrogen)atoms.getAtom(2)).getOrientation().getDirection()));
//            if (Math.abs(phi02new - newAlpha[2]) > 1.0E-10 ) {
//                System.out.println(phi02new+" "+newAlpha[2]);
//                throw new RuntimeException();
//            }
//            if (printToFile) {
//                try {
//                    newFancy.write(phi02new+" ");
//                }
//                catch (IOException e) {
//                    throw new RuntimeException(e);
//                }            
//            }
//            break;
        }
        double pActRatio = Math.exp(uOld-uNew);
        uAcc = uNew;
        pacc = pActRatio/pGenRatio;
        
        double pActOld = Math.exp(-uOld);
        double pActNew = Math.exp(-uNew);
//        System.out.println(pacc+" "+pGenOld+" "+pActOld+" "+pGenNew+" "+pActNew);
//        System.out.println(pGenNew);
//        if (printToFile) {
//            try {
//                oldFancy.write(pGenOld+" "+pActOld+"\n");
//                newFancy.write(pGenNew+" "+pActNew+"\n");
//            }
//            catch (IOException e) {
//                throw new RuntimeException(e);
//            }
//        }
//        if (box.getIndex() == 0) {
//            System.out.println("fancy: "+uOld+" "+uNew);        
//            if (uCurrent == 0) uCurrent = uOld;

//            if (Math.abs(uOld - uCurrent) > 1.0E-10) throw new RuntimeException();
//        }
//        System.out.println("pacc = "+pacc);
//        flag = false;
		((BoxCluster)box).trialNotify();
        weightNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
       
		return true;
	}
	protected FileWriter oldFancy,newFancy;
    public boolean printToFile;
    public double getA() {
        if (firstMove) {
            return weightNew/weightOld;
        }
        else {
//            System.out.println("pacc = "+pacc);
//            if (weightNew != weightOld) throw new RuntimeException("busted");
            return pacc*weightNew/weightOld;
        }
//        we skip calculation of wOld because we're the only intramolecular move in town.
//        System.out.println("return = "+pacc*weightNew/weightOld);
//        return 1.00;
    }

    public double getB() {
    	return 0.0;
    }
    
    public void rejectNotify() {
        rejected++;
        prevAcc = false;
//        if (!flag) hist1.addValue(newAlpha[2], 0.0);            
//        System.out.println("rejected");
        IMoleculeList molecules = box.getMoleculeList();
        for (int i=0; i<molecules.getMoleculeCount(); i++) {
            IAtomList atoms = molecules.getMolecule(i).getChildList();            
            for (int k=0; k<P; k++) {                
                    IAtomOriented kAtom = ((IAtomOriented)atoms.getAtom(k));                    
                    kAtom.getOrientation().setDirection(oldOrientations[i][k].getDirection());                    
            }
        }
    	((BoxCluster)box).rejectNotify();
    }
    
    public void acceptNotify() {
//    	System.out.println("orientation move accepted");
        if (box.getIndex() == 0) uCurrent = uAcc;
        if (rejected > 100 && box.getIndex() == 0) {
//            System.out.println("Rejected "+rejected+" moves continuously, printing to Orientation"+rc+".dat");
//            try {
//                file = new FileWriter("Orientation"+rc+".dat");
//                for (int i = 0; i<box.getMoleculeList().getMoleculeCount(); i++) {
//                    for (int j = 0; j<P; j++) {
//                        file.write(((AtomHydrogen)box.getMoleculeList().getMolecule(i).getChildList().getAtom(j)).getBondLength()+", "+oldOrientations[i][j].getDirection()+"\n");
//                    }
//                }
//                file.close();
//            }
//            catch (IOException e) {
//                throw new RuntimeException(e);
//            }
            rc++;
//            System.exit(1);
        }
        prevAcc = true;
        rejected = 0;
        foo++;
        firstMove = false;
//        hist1.addValue(newAlpha[2], 1);
//        hist2.addValue(newAlpha[2]);
//        flag = true;        
//        System.out.println("fancy move accepted");
//    	if (firstMove) {
//    	    ((BoxCluster)box).acceptNotify();
//    	    firstMove = false;
//    	}
//    	else {
//    	    rejectNotify();
//    	}
        IMoleculeList molecules = box.getMoleculeList();
        for (int i=0; i<molecules.getMoleculeCount(); i++) {                        
            for (int k=0; k<P; k++) {
                xOld[i][k] = xNew[i][k];
            }
        }            
        ((BoxCluster)box).acceptNotify();
    }
    
    public double energyChange() {
        return 0;
    }
    
    public AtomIterator affectedAtoms() {
        return leafIterator;
    }
    protected boolean doExchange = false;
    public void setDoExchange(boolean b) {
    	doExchange = b;
    }
    protected double dist(AtomHydrogen a0, AtomHydrogen a1) {
    	dummy.E(a0.getOrientation().getDirection());
    	if (a0.getIndex() == 0 && doExchange) {
    		dummy.TE(-1);
    	}
        dr.Ev1Mv2(a0.getPosition(),a1.getPosition());
        dr1.Ea1Tv1(a0.getBondLength()/2, dummy);
        dr1.PEa1Tv1(-a1.getBondLength()/2, a1.getOrientation().getDirection());
//        if (a0.getIndex() == 0 && doExchange) a0.getOrientation().setDirection(dummy);
        double r2 = dr.Mv1Squared(dr1);
//        if (flag) {
//            System.out.print("r = "+Math.sqrt(r2)+" ");
//            flag = false;
//        }
        dr1.TE(-1);
        r2 += dr.Mv1Squared(dr1);
        return r2;
    }
    public static void rotateBy(double cdt, double sdt, IVector axis, IVectorMutable direction) {
        // consider a circle on the surface of the unit sphere.  The given axis
        // passes through the center of the circle.  The circle passes through
        // the current direction vector and the vector v4 defined below.  We
        // rotate the direction by the given angle (dt) around the circle.
        
        // v1 is the projection of direction onto axis
        // v2 is the component of direction perpendicular to axis
        // v3 has the same magnitude as v2 and is perpendicular to both
        //    direction and axis
        // v4 is a unit vector whose components are v1 and v3
        
        // v1 = v1overAxis * axis
        ISpace space = Space3D.getInstance();        
        double v1overAxis = axis.dot(direction);
        IVectorMutable temp = space.makeVector();
        IVectorMutable temp2 = space.makeVector();
        temp.Ea1Tv1(-v1overAxis, axis);
        temp.PE(direction);
        // now temp = v2
        temp2.E(axis);
        temp2.XE(direction);
        // now temp2 = v3
        direction.Ea1Tv1(cdt, temp);
        direction.PEa1Tv1(sdt, temp2);
        direction.PEa1Tv1(v1overAxis, axis);
        
    }
    protected double uAcc = 0;
    protected FileWriter file;
    public static double uCurrent = 0;
    private static final long serialVersionUID = 1L;
    protected final int P;    
    protected final ISpace space;
    protected final IRandom random;    
    protected IOrientation3D[][] oldOrientations;    
    // cluster weights
    protected double weightOld, weightNew;    
    protected final AtomIteratorLeafAtoms leafIterator;
    protected int foo = 0;
    protected double stiffness;
    protected double pacc;
    protected boolean prevAcc = false;
    protected boolean firstMove = true;
    protected int rejected = 0;
    protected int rc = 0;
    protected int choice = 1;
    public HistogramNotSoSimple hist1 = new HistogramNotSoSimple(10000,new DoubleRange(0, 1));
    public HistogramSimple hist2 = new HistogramSimple(10000,new DoubleRange(0, 1));
    protected double [] newAlpha;
    protected double [][] xOld;
    protected double [][] xNew;
    protected boolean flag;
    protected final IVectorMutable dr,dr1,dr2,dummy;
    public static void main(String[] args) {
        ISpace space = Space3D.getInstance();
        ClusterWeight cluster = new ClusterWeight() {
            
            @Override
            public double value(BoxCluster box) {
                // TODO Auto-generated method stub
                return 1;
            }
            
            @Override
            public void setTemperature(double temperature) {
                // TODO Auto-generated method stub
                
            }
            
            @Override
            public int pointCount() {
                // TODO Auto-generated method stub
                return 1;
            }
            
            @Override
            public ClusterAbstract makeCopy() {
                // TODO Auto-generated method stub
                return null;
            }
        };
        BoxCluster box = new BoxCluster(cluster,space);
        Simulation sim = new Simulation(space);
        sim.addBox(box);
        IAtomTypeOriented atype = new AtomTypeOrientedSphere(Hydrogen.INSTANCE, space);
        SpeciesSpheresHetero species = new SpeciesSpheresHetero(space,new IAtomTypeOriented [] {atype});
        sim.addSpecies(species);
        File file1 = new File ("acceptance.dat");
        if (file1.exists()) {
            file1.delete();
        }
        for (int p = 2; p <= 512; p*=2) {
            box.setNMolecules(species,0);
            species.setChildCount(new int [] {p});
            box.setNMolecules(species, 1);
            IntegratorMC integrator = new IntegratorMC(sim, null);
            integrator.setBox(box);
            MCMoveClusterRingRegrowOrientation move = new MCMoveClusterRingRegrowOrientation(sim.getRandom(), space, p);
            
            for (int iTemp = 40; iTemp <= 40; iTemp+= 2) {               
                move.foo = 0;
                move.setStiffness(Kelvin.UNIT.toSim(iTemp), species.getAtomType(0).getMass());
                integrator.getMoveManager().addMCMove(move);
                integrator.reset();
                int total = 100;
                for (int i=0; i<total; i++) {
                    integrator.doStep();
                }
                        try{                        
                            FileWriter Temp = new FileWriter("acceptance.dat",true);
                            Temp.write(iTemp+" "+p+" "+move.getStiffness()+" "+((double)move.foo)/total+"\n");                          
                            Temp.close();
                        }
                        catch(IOException ex1){
                            throw new RuntimeException(ex1);
                        }
                System.out.println("p = "+p+" ,Temp = "+iTemp+" ,acceptance ratio = "+((double)move.foo)/total);
            }
            
        }
    }
}


