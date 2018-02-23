/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.action.MoleculeActionTranslateTo;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.box.BoxAgentManager;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorRigidIterative;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisHcp;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveHexagonal;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculePositionDefinition;
import etomica.molecule.MoleculeAgentManager;
import etomica.molecule.MoleculeAgentManager.MoleculeAgentSource;
import etomica.molecule.MoleculePositionGeometricCenter;
import etomica.nbr.cell.molecule.NeighborCellManagerMolecular;
import etomica.nbr.list.molecule.BoxAgentSourceCellManagerListMolecular;
import etomica.nbr.list.molecule.NeighborListManagerSlantyMolecular;
import etomica.nbr.list.molecule.PotentialMasterListMolecular;
import etomica.normalmode.BasisBigCell;
import etomica.normalmode.MeterHarmonicEnergy;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationEnergySumBigDecimal;
import etomica.potential.PotentialCalculationTorqueSum;
import etomica.simulation.Simulation;
import etomica.space.*;
import etomica.space3d.RotationTensor3D;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.units.Degree;

import java.util.Arrays;



/**
 * Lattice-Energy Minimization for Beta-Nitrogen Structure
 * using gradient
 *  
 * @author Andrew Schultz
 *
 */
public class MinimizationBetaNitrogenModel extends Simulation{

	public MinimizationBetaNitrogenModel(Space space, int[] nC, double density) {
        super(space);
        this.space = space;

        BoxAgentSourceCellManagerListMolecular boxAgentSource = new BoxAgentSourceCellManagerListMolecular(this, null, space);
        BoxAgentManager<NeighborCellManagerMolecular> boxAgentManager = new BoxAgentManager<NeighborCellManagerMolecular>(boxAgentSource, this);

        double ratio = 1.631;
        double a = Math.pow(4.0 / (Math.sqrt(3.0) * ratio * density), 1.0 / 3.0);
        double c = a * ratio;
        int numMolecule = nC[0] * nC[1] * nC[2] * 2;
        int nCa = (int) Math.pow(numMolecule / 1.999999999, 1.0 / 3.0);

        Basis basisHCP = new BasisHcp();
        Basis basis = new BasisBigCell(space, basisHCP, new int[]{nC[0], nC[1], nC[2]});

        species = new SpeciesN2(space);
        addSpecies(species);

        box = this.makeBox();
        box.setNMolecules(species, numMolecule);
        int[] nCells = new int[]{1, 1, 1};

        Vector[] boxDim = new Vector[3];
        boxDim[0] = space.makeVector(new double[]{nC[0] * a, 0, 0});
        boxDim[1] = space.makeVector(new double[]{-nC[1] * a * Math.cos(Degree.UNIT.toSim(60)), nC[1] * a * Math.sin(Degree.UNIT.toSim(60)), 0});
        boxDim[2] = space.makeVector(new double[]{0, 0, nC[2] * c});

        Boundary boundary = new BoundaryDeformablePeriodic(space, boxDim);
        primitive = new PrimitiveHexagonal(space, (nC[0]) * a, nC[2] * c);

        coordinateDef = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
        coordinateDef.setIsBeta();
        coordinateDef.setOrientationVectorBeta(space);
        coordinateDef.initializeCoordinates(nCells);


        double[] u = new double[20];

        if (true) {
            BetaPhaseLatticeParameter parameters = new BetaPhaseLatticeParameter();
            double[][] param = parameters.getParameter(density);

            int kParam = 0;
            for (int i = 0; i < param.length; i++) {
                for (int j = 0; j < param[0].length; j++) {
                    u[kParam] = param[i][j];
                    kParam++;
                }
            }

            int numDOF = coordinateDef.getCoordinateDim();
            double[] newU = new double[numDOF];

//			double[] deviation = new double[]{
//					1.2780887459484802E-10, -7.492459985769528E-10, 1.4581758023268776E-10, 2.4394273825234006E-8, 8.42930716693827E-9, 
//					-1.3034284762625248E-10, 7.715463823387836E-10, -1.397282289872237E-10, -2.4096713784662293E-8, -9.245658416530434E-9, 
//					-1.4337864229219122E-10, -7.749569874704321E-10, -1.4337686593535182E-10, 3.1120726966510594E-8, -1.1904762815825307E-8, 
//					1.4595613606616098E-10, 7.526033130034193E-10, 1.3732659454035456E-10, -0.0, 0.0
//			};
//			
//			for(int i=0; i<u.length; i++){
//				u[i] += deviation[i];
//				System.out.print(u[i]+", ");
//				if(i%5==4) System.out.println();
//			}
//			
//			System.exit(1);
            if (true) {
                for (int j = 0; j < numDOF; j += 10) {
                    if (j > 0 && j % (nCa * 10) == 0) {
                        j += nCa * 10;
                        if (j >= numDOF) {
                            break;
                        }
                    }
                    for (int k = 0; k < 10; k++) {
                        newU[j + k] = u[k];
                    }
                }

                for (int j = nCa * 10; j < numDOF; j += 10) {
                    if (j > nCa * 10 && j % (nCa * 10) == 0) {
                        j += nCa * 10;
                        if (j >= numDOF) {
                            break;
                        }
                    }
                    for (int k = 0; k < 10; k++) {
                        newU[j + k] = u[k + 10];
                    }
                }
            }

            coordinateDef.setToU(box.getMoleculeList(), newU);
            coordinateDef.initNominalU(box.getMoleculeList());

        }
        this.initialU = u;

        box.setBoundary(boundary);
        double rC = a * nC[0] * 0.475;
        //System.out.println("Truncation Radius: " + rC);
        potential = new P2Nitrogen(space, rC);
        potential.setBox(box);

        potentialMaster = new PotentialMasterListMolecular(this, rC, boxAgentSource, boxAgentManager, new NeighborListManagerSlantyMolecular.NeighborListSlantyAgentSourceMolecular(rC, space), space);
        potentialMaster.addPotential(potential, new ISpecies[]{species, species});
        int cellRange = 6;
        potentialMaster.setRange(rC);
        potentialMaster.setCellRange(cellRange);
        potentialMaster.getNeighborManager(box).reset();

        potential.setRange(Double.POSITIVE_INFINITY);
        int potentialCells = potentialMaster.getNbrCellManager(box).getLattice().getSize()[0];
        if (potentialCells < cellRange * 2 + 1) {
            throw new RuntimeException("oops (" + potentialCells + " < " + (cellRange * 2 + 1) + ")");
        }

    }
	
	public static void main (String[] args){

		int nC0 = 8; 
		int nC1 = 8; 
		int nC2 = 8;
		int[] nC = new int []{nC0,nC1,nC2};
		int numMolecule =nC[0]*nC[1]*nC[2]*2;

		double density = 0.0210;
		final MinimizationBetaNitrogenModel sim = new MinimizationBetaNitrogenModel(Space3D.getInstance(3), nC, density);
	    
		final MeterPotentialEnergy meterPotentialEnergy = new MeterPotentialEnergy(sim.potentialMaster, sim.box);
		meterPotentialEnergy.setPotentialCalculation(new PotentialCalculationEnergySumBigDecimal(18));
	
		double latticeEnergy = meterPotentialEnergy.getDataAsScalar();
	
		PotentialCalculationTorqueSum pcForce = new PotentialCalculationTorqueSum();
        MoleculeAgentSource molAgentSource = new MoleculeAgentSource() {
            
            public void releaseAgent(Object agent, IMolecule molecule) {
            }
            
            public Object makeAgent(IMolecule mol) {
                return new IntegratorRigidIterative.MoleculeAgent(sim.space);
            }

        };
        
        MoleculeAgentManager molAgentManager = new MoleculeAgentManager(sim, sim.box, molAgentSource);
        pcForce.setMoleculeAgentManager(molAgentManager);
        double[] d = new double[16];
        MoleculeActionTranslateTo translator = new MoleculeActionTranslateTo(sim.space);
        MoleculePositionGeometricCenter pos = new MoleculePositionGeometricCenter(sim.space);
        translator.setAtomPositionDefinition(pos);
        Vector p = sim.space.makeVector();
        int nA = 16;
        double step1 = 0;
        double[] x0 = new double[12];
        
        Vector[] orient0 = new Vector[4];
        Vector[] orientf = new Vector[4];
        Vector[] torques = new Vector[4];
        Vector[] axes = new Vector[4];
        for (int i=0; i<4; i++) {
            axes[i] = sim.space.makeVector();
            torques[i] = sim.space.makeVector();
            orient0[i] = sim.space.makeVector();
            orientf[i] = sim.space.makeVector();
        }
        for (int outer = 0; outer < 20; outer++) {
            System.out.println("**** "+outer+" ****");
	        double totalD = 0;
	        double step = 0;
	        double radianFac = 1; //????

	        RotationTensor3D rTensor = new RotationTensor3D();

		    for (int iter=0; iter<3; iter++) {
		        double[] g = new double[12];
		        double t = 0;
		        for (int i=0; i<4; i++) {
		            IMolecule iMol = sim.box.getMoleculeList().get(i<2 ? i : i+(nA-2));
		            IteratorDirective id = new IteratorDirective(null, iMol);
		            sim.potentialMaster.calculate(sim.box, id, pcForce);
		            Vector f = ((IntegratorRigidIterative.MoleculeAgent)molAgentManager.getAgent(iMol)).force;
		            for (int j=0; j<3; j++) {
		                g[i*3+j] = -f.getX(j);
		                t += g[i*3+j]*g[i*3+j];
		            }
		            torques[i].E(((IntegratorRigidIterative.MoleculeAgent)molAgentManager.getAgent(iMol)).torque);
//		            System.out.println(torques[i]);
		            
		            if (iter == 0) {
		                double t2 = torques[i].squared();
		                double sqrtT = Math.sqrt(t2);
		                t += t2/(radianFac*radianFac);
		                d[12+i] = sqrtT/radianFac;
		                axes[i].E(torques[i]);
		                axes[i].TE(1/sqrtT);
		            }
		            pcForce.reset();
		        }
		        System.out.println(Arrays.toString(g));
		        if (iter == 0) {
		            t = Math.sqrt(t);

		            for (int j=0; j<12; j++) {
		                d[j] = -g[j]/t;
		                totalD += g[j]*d[j];
		            }
		            for (int j=12; j<16; j++) {
		                d[j] /= t;
		                totalD -= d[j]*torques[j-12].dot(axes[j-12])/radianFac;
		            }
		            System.out.println(Arrays.toString(d));
		            System.out.println("totalD "+totalD);
		        
		            if (step1 == 0) {
		                step = 0.0000001;
		            }
		            else {
		                step = step1/5;
		            }
		            System.out.println("step "+step);
		        }
		        else {
		        	double newTotalD = 0;
	                for (int j=0; j<12; j++) {
	                    newTotalD += g[j]*d[j];
	                }

	                for (int j=12; j<16; j++) {
	                	newTotalD -= d[j]*torques[j-12].dot(axes[j-12])/radianFac;
	                }
	                
	                System.out.println("totalD "+newTotalD);
	                

	                double oldStep = step; 
	                if(newTotalD > totalD || totalD > 0.0){
	                	
	                	step = - newTotalD * step / (newTotalD - totalD);
	                	if (iter == 1) {
		                    step1 = oldStep + step;
		                }
	                } else {
	                	
	                	step = step *2;
	                	if (iter == 1) {
		                    step1 = step*4;
		                }
	                }

	                totalD = newTotalD;
	                System.out.println("step "+step);
		        }
		
		        for (int i=0; i<4; i++) {
		            for (int k=i%2; k<numMolecule; k+=2) {
		                boolean isA = (k/nA)%2 == 0;
		                if ((i<2 && !isA) || (i>1 && isA)) continue;
		//                    System.out.println(i+" "+k);
		                IMolecule iMol = sim.box.getMoleculeList().get(k);
		                p.E(pos.position(iMol));
		                for (int j=0; j<3; j++) {
		                    if (x0[i*3+j] == 0) {
		                        x0[i*3+j] = p.getX(j);
		                    }
		                    p.setX(j, p.getX(j)+step*d[i*3+j]);
		                }
		                translator.setDestination(p);
		                translator.actionPerformed(iMol);
		                
		                if (orient0[i].isZero()) {
	                        orient0[i].E(iMol.getChildList().get(0).getPosition());
	                        orient0[i].ME(pos.position(iMol));
	                        orient0[i].normalize();
	                    }
		                
		                rTensor.setRotationAxis(axes[i], step*d[12+i]/radianFac);
		                doTransform(iMol, pos, rTensor);
		            }
		        }
		    }
        }
        

        //DONE with minimization
        
        double[] newU = sim.coordinateDef.calcU(sim.box.getMoleculeList());
        double[] deviation = new double[20];
        
        System.out.println("\ndeviation in newU");
        for(int i=0; i<4; i++){
        	int iMolec = i<2 ? i : i+(nA-2);
        	for (int k=0; k<5; k++){
        		deviation[i*5+k] = newU[iMolec*5+k];
        		System.out.print(newU[iMolec*5+k]+", ");
        	}
        	System.out.println();
        }
    
        double[] xf = new double[12];
        double disp = 0.0;
        double angleDisp = 0.0;
        for (int i=0; i<4; i++) {
            IMolecule iMol = sim.box.getMoleculeList().get(i<2 ? i : i+(nA-2));
            p.E(pos.position(iMol));
            for (int j=0; j<3; j++) {
//                System.out.println(x0[i*3+j]+" => "+p.getX(j)+"    "+(p.getX(j)-x0[i*3+j]));
                xf[i*3+j] = p.getX(j);
                double dx = xf[i*3+j] - x0[i*3+j];
                disp += dx*dx;
            }
            
            orientf[i].E(iMol.getChildList().get(0).getPosition());
            orientf[i].ME(pos.position(iMol));
            orientf[i].normalize();
            angleDisp += orientf[i].Mv1Squared(orient0[i]);
//            System.out.println("     "+Math.sqrt(orientf[i].Mv1Squared(orient0[i])));
        }
        
        disp = Math.sqrt(disp);
        angleDisp = Math.sqrt(angleDisp);
        System.out.println("\ndisp "+disp+"  angleDisp "+angleDisp);
        double newLatticeEnergy = meterPotentialEnergy.getDataAsScalar();
        System.out.println("Old Lattice Energy (per molecule): "+latticeEnergy/numMolecule);
        System.out.println("New Lattice Energy (per molecule): "+newLatticeEnergy/numMolecule);
        
		for(int i=0; i<initialU.length; i++){
			if(i%5==0) System.out.print("{");
			initialU[i] += deviation[i];
			System.out.print(initialU[i]);
			if(i%5==4) {
				System.out.println("},");
			} else {
				System.out.print(", ");
			}
		}
        System.exit(1);
        for (int l=0; l<201; l++) {
            for (int i=0; i<4; i++) {
                for (int k=i%2; k<numMolecule; k+=2) {
                    boolean isA = (k/nA)%2 == 0;
                    if ((i<2 && !isA) || (i>1 && isA)) continue;
                    IMolecule iMol = sim.box.getMoleculeList().get(k);
                    p.E(pos.position(iMol));
                    for (int j=0; j<3; j++) {
                        if (l==0) {
                            p.setX(j, p.getX(j)-(xf[i*3+j] - x0[i*3+j])*10);
                        }
                        else {
                            p.setX(j, p.getX(j)+(xf[i*3+j] - x0[i*3+j])/10.0);
                        }
                    }
                    translator.setDestination(p);
                    translator.actionPerformed(iMol);
                }
            }
            double u = 0;
            for (int i=0; i<4; i++) {
                for (int k=(4*nA)+i%2; k<numMolecule; k+=2) {
                    boolean isA = (k/nA)%2 == 0;
                    if ((i<2 && !isA) || (i>1 && isA)) continue;
                    IMolecule iMol = sim.box.getMoleculeList().get(k);
                    meterPotentialEnergy.setTarget(iMol);
                    u += meterPotentialEnergy.getDataAsScalar();
                    break;
                }
            }
//            meterPotentialEnergy.setTarget((IMolecule)null);
//            System.out.println((l-100)/100.0+" "+(meterPotentialEnergy.getDataAsScalar()-newLatticeEnergy)/numMolecule);
            System.out.println((l-100)/100.0+" "+(u/8-newLatticeEnergy/numMolecule));
        }            
        System.exit(1);

	}
	
    protected static void doTransform(IMolecule molecule, IMoleculePositionDefinition posDef, Tensor rotationTensor) {
        IAtomList childList = molecule.getChildList();
        Vector com = posDef.position(molecule);
        for (int iChild = 0; iChild<childList.size(); iChild++) {
            IAtom a = childList.get(iChild);
            Vector r = a.getPosition();
            r.ME(com);
            rotationTensor.transform(r);
            r.PE(com);
        }
    }


	protected Box box;
	protected Space space;
	protected PotentialMasterListMolecular potentialMaster;
	protected P2Nitrogen potential;
	protected CoordinateDefinitionNitrogen coordinateDef;
	protected Primitive primitive;
	protected SpeciesN2 species;
	protected MeterHarmonicEnergy meterHarm;
	protected MeterPotentialEnergy meterPE;
	protected static double[] initialU;
	private static final long serialVersionUID = 1L;
}
