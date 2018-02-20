/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.action.AtomActionTranslateBy;
import etomica.action.MoleculeActionTranslateTo;
import etomica.action.MoleculeChildAtomAction;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataInfo;
import etomica.data.FunctionData;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisHcp;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveTriclinic;
import etomica.models.nitrogen.LatticeSumCrystalMolecular.DataGroupLSC;
import etomica.molecule.*;
import etomica.normalmode.BasisBigCell;
import etomica.paracetamol.AtomActionTransformed;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space3d.RotationTensor3D;
import etomica.space3d.Space3D;
import etomica.units.Degree;
import etomica.units.dimensions.Energy;

import java.util.Arrays;



/**
 * Lattice-Energy Minimization for Beta-Nitrogen Structure
 * using gradient with lattice sum
 *  
 * @author Andrew Schultz and Tai Boon Tan
 *
 */
public class MinimizationBetaNitrogenModelLS extends Simulation{

	public MinimizationBetaNitrogenModelLS(Space space, double density, double rC) {
        super(space);
        this.space = space;
        this.density = density;

        double ratio = 1.631;
        double aDim = Math.pow(4.0 / (Math.sqrt(3.0) * ratio * density), 1.0 / 3.0);
        double cDim = aDim * ratio;
        System.out.println("density: " + density);
        System.out.println("aDim: " + aDim + " ;cDim: " + cDim);

        int[] nCells = new int[]{1, 2, 1};
        Basis basisHCP = new BasisHcp();
        basis = new BasisBigCell(space, basisHCP, nCells);

        ConformationNitrogen conformation = new ConformationNitrogen(space);
        SpeciesN2 species = new SpeciesN2(space);
        species.setConformation(conformation);
        addSpecies(species);

        SpeciesN2B ghostSpecies = new SpeciesN2B(space);
        ghostSpecies.setConformation(conformation);
        addSpecies(ghostSpecies);

        int numMolecule = 4;
        box = this.makeBox();
        box.setNMolecules(species, numMolecule);

        ghostBox = this.makeBox();
        ghostBox.setNMolecules(ghostSpecies, 1);

        Vector[] boxDim = new Vector[3];
        boxDim[0] = space.makeVector(new double[]{aDim, 0, 0});
        boxDim[1] = space.makeVector(new double[]{-2 * aDim * Math.cos(Degree.UNIT.toSim(60)), 2 * aDim * Math.sin(Degree.UNIT.toSim(60)), 0});
        boxDim[2] = space.makeVector(new double[]{0, 0, cDim});

        BoundaryDeformablePeriodicSwitch boundary = new BoundaryDeformablePeriodicSwitch(space, boxDim);
        boundary.setDoPBC(false);
        primitive = new PrimitiveTriclinic(space, aDim, 2 * aDim, cDim, Math.PI * (90 / 180.0), Math.PI * (90 / 180.0), Math.PI * (120 / 180.0));
        box.setBoundary(boundary);

        BetaPhaseLatticeParameterLS parameters = new BetaPhaseLatticeParameterLS();
        double[][] param = parameters.getParameter(density);

        double[] u = new double[20];
        int kParam = 0;
        for (int i = 0; i < param.length; i++) {
            for (int j = 0; j < param[0].length; j++) {
                u[kParam] = param[i][j];
                kParam++;
            }
        }
        this.u = u;

        coordinateDef = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
        coordinateDef.setIsBetaLatticeSum();
        coordinateDef.setIsDoLatticeSum();
        coordinateDef.setOrientationVectorBetaLatticeSum(space, density, param);
        coordinateDef.initializeCoordinates(new int[]{1, 1, 1});

        potential = new P2Nitrogen(space, rC);
        potential.setBox(box);
        potential.setEnablePBC(false);

        this.nLayer = (int) (rC / aDim + 0.5);

        FunctionData<Object> function = new FunctionData<Object>() {
            public IData f(Object obj) {
                data.x = potential.energy((IMoleculeList) obj);
                return data;
            }

            public IDataInfo getDataInfo() {
                return dataInfo;
            }

            final DataInfo dataInfo = new DataDouble.DataInfoDouble("Lattice energy", Energy.DIMENSION);
            final DataDouble data = new DataDouble();
        };

        BravaisLatticeCrystal lattice = new BravaisLatticeCrystal(primitive, basis);
        LatticeSumCrystalMolecular latticeSum = new LatticeSumCrystalMolecular(lattice, coordinateDef, ghostBox);
        latticeSum.setMaxLatticeShell(nLayer);

        double sum = 0;
        double basisDim = lattice.getBasis().getScaledCoordinates().length;
        DataGroupLSC data = (DataGroupLSC) latticeSum.calculateSum(function);
        for (int j = 0; j < basisDim; j++) {
            for (int jp = 0; jp < basisDim; jp++) {
                sum += ((DataDouble) data.getDataReal(j, jp)).x;
            }
        }

        initialLatticeEnergy = (0.5 * sum / basisDim);
        System.out.println("initial energy: " + initialLatticeEnergy);


        atomGroupAction = new MoleculeChildAtomAction(new AtomActionTransformed(lattice.getSpace()));

        translateBy = new AtomActionTranslateBy(coordinateDef.getPrimitive().getSpace());
        atomGroupActionTranslate = new MoleculeChildAtomAction(translateBy);
        lsPosition = space.makeVector();

        xVecBox = Math.sqrt(box.getBoundary().getEdgeVector(0).squared());
        yVecBox = Math.sqrt(box.getBoundary().getEdgeVector(1).squared());
        zVecBox = Math.sqrt(box.getBoundary().getEdgeVector(2).squared());

    }
	
	public double getEnergy (double[] u){

		double param[][] = new double[4][5];
		for(int i=0; i<4; i++){
			for(int j=0; j<5; j++){
				param[i][j] = u[i*5+j];
					
			}	
		}
		
		CoordinateDefinitionNitrogen coordDef = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
		coordDef.setIsBetaLatticeSum();
		coordDef.setIsDoLatticeSum();
		coordDef.setOrientationVectorBetaLatticeSum(space, density, param);
		coordDef.initializeCoordinates(new int[]{1,1,1});
		
		FunctionData<Object> function = new FunctionData<Object>() {
			public IData f(Object obj) {
				data.x = potential.energy((IMoleculeList)obj);
				return data;
			}
			public IDataInfo getDataInfo() {
				return dataInfo;
			}
			final DataInfo dataInfo = new DataDouble.DataInfoDouble("Lattice energy", Energy.DIMENSION);
			final DataDouble data = new DataDouble();
		};
		
		BravaisLatticeCrystal lattice = new BravaisLatticeCrystal(primitive, basis);
		LatticeSumCrystalMolecular latticeSum = new LatticeSumCrystalMolecular(lattice, coordDef, ghostBox);
		latticeSum.setMaxLatticeShell(nLayer);
		
		double sum = 0;
	    double basisDim = lattice.getBasis().getScaledCoordinates().length;
		DataGroupLSC data = (DataGroupLSC)latticeSum.calculateSum(function);
        for(int j=0; j<basisDim; j++) {
            for(int jp=0; jp<basisDim; jp++) {
                sum += ((DataDouble)data.getDataReal(j,jp)).x; 
            }
        }

		return 0.5*sum/basisDim;
	}
	
	public static void main (String[] args){
	
		double rC = 250.0;
		double density = 0.0223;
		if(args.length > 0){
			density = Double.parseDouble(args[0]);
		}
		
		final MinimizationBetaNitrogenModelLS sim = new MinimizationBetaNitrogenModelLS(Space3D.getInstance(3), density, rC);
	
		
        double[] d = new double[16];
        MoleculeActionTranslateTo translator = new MoleculeActionTranslateTo(sim.space);
        MoleculePositionGeometricCenter pos = new MoleculePositionGeometricCenter(sim.space);
        translator.setAtomPositionDefinition(pos);
        Vector p = sim.space.makeVector();

        double step1 = 0;
        double[] x0 = new double[12];
        
        Vector[][] gradTorq = new Vector[2][2];
        for(int i=0; i<2; i++){
        	for(int j=0; j<2; j++){
            	gradTorq[i][j] = sim.getSpace().makeVector();
            }	
        }
        
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
		        
		        //Looping over molecules
	            //To get the gradient from force calculation
	            MoleculePair pair = new MoleculePair();
		        	            
		        for (int i=0; i<4; i++) {
		            IMolecule molecule0 = sim.box.getMoleculeList().getMolecule(i);
		            pair.atom0 = molecule0;
		            
		            Vector destination = sim.getSpace().makeVector();
		            Vector f = sim.getSpace().makeVector();
		            torques[i].E(0.0);
		            
		           	for(int jp=0; jp<4; jp++){	
		                IMolecule ghostMol = sim.ghostBox.getMoleculeList(sim.ghostBox.getMoleculeList().getMolecule(0).getType()).getMolecule(0);
	                    ghostMol.getType().initializeConformation(ghostMol);
	                    
	                    pair.atom1 = ghostMol;
	            		
	                    int rotationNum = jp%4;
	                    ((AtomActionTransformed)sim.atomGroupAction.getAtomAction()).setTransformationTensor(sim.coordinateDef.yOrientationTensor[rotationNum]);
		                sim.atomGroupAction.actionPerformed(ghostMol);
		                
		                ((AtomActionTransformed)sim.atomGroupAction.getAtomAction()).setTransformationTensor(sim.coordinateDef.xzOrientationTensor[rotationNum]);
		                sim.atomGroupAction.actionPerformed(ghostMol);
		            	
		            	destination.E(pos.position(sim.box.getMoleculeList().getMolecule(jp)));
		        		translator.setDestination(destination);
						translator.actionPerformed(pair.atom1); 
						
		            	int nLayer = sim.nLayer;
		            	for(int x=-nLayer; x<=nLayer; x++){
		    				for(int y=-nLayer; y<=nLayer; y++){
		    					for(int z=-nLayer; z<=nLayer; z++){
		    						if(x==0 && y==0 && z==0 && (i==jp)) continue;
		    						sim.lsPosition.E(new double[]{x*sim.xVecBox-y*sim.yVecBox*Math.cos(Degree.UNIT.toSim(60)), 
		    													 y*sim.yVecBox*Math.sin(Degree.UNIT.toSim(60)), 
		    													 z*sim.zVecBox});

		    						sim.translateBy.setTranslationVector(sim.lsPosition);
		    						sim.atomGroupActionTranslate.actionPerformed(ghostMol);
		    						
		    						gradTorq = sim.potential.gradientAndTorque(pair);
		    						f.ME(gradTorq[0][0]);
		    						torques[i].PE(gradTorq[1][0]);
		    						
		    						sim.lsPosition.TE(-1);
		    						sim.translateBy.setTranslationVector(sim.lsPosition);
		    						sim.atomGroupActionTranslate.actionPerformed(ghostMol);
		    					}	
		    				}	
		    			}
		           }
		            
		            for (int j=0; j<3; j++) {
		                g[i*3+j] = -f.getX(j);
		                t += g[i*3+j]*g[i*3+j];
		            }
//		            torques[i].E(((IntegratorRigidIterative.MoleculeAgent)molAgentManager.getAgent(iMol)).torque);
//		            System.out.println(torques[i]);
		            
		            if (iter == 0) {
		                double t2 = torques[i].squared();
		                double sqrtT = Math.sqrt(t2);
		                t += t2/(radianFac*radianFac);
		                d[12+i] = sqrtT/radianFac;
		                axes[i].E(torques[i]);
		                axes[i].TE(1/sqrtT);
		            }
		            
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
		        } else {
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
		            IMolecule iMol = sim.box.getMoleculeList().getMolecule(i);
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
                        orient0[i].E(iMol.getChildList().getAtom(0).getPosition());
                        orient0[i].ME(pos.position(iMol));
                        orient0[i].normalize();
                    }
	                
	                rTensor.setRotationAxis(axes[i], step*d[12+i]/radianFac);
	                doTransform(iMol, pos, rTensor);
		        }
		    }
        }
        

        //DONE with minimization
        
        double[] newU = sim.coordinateDef.calcU(sim.box.getMoleculeList());
        double[] deviation = new double[20];
        
        System.out.println("\ndeviation in newU");
        for(int i=0; i<4; i++){
        	for (int k=0; k<5; k++){
        		deviation[i*5+k] = newU[i*5+k];
        		System.out.print(newU[i*5+k]+", ");
        	}
        	System.out.println();
        }
    
        double[] xf = new double[12];
        double disp = 0.0;
        double angleDisp = 0.0;
        for (int i=0; i<4; i++) {
            IMolecule iMol = sim.box.getMoleculeList().getMolecule(i);
            p.E(pos.position(iMol));
            for (int j=0; j<3; j++) {
//                System.out.println(x0[i*3+j]+" => "+p.getX(j)+"    "+(p.getX(j)-x0[i*3+j]));
                xf[i*3+j] = p.getX(j);
                double dx = xf[i*3+j] - x0[i*3+j];
                disp += dx*dx;
            }
            
            orientf[i].E(iMol.getChildList().getAtom(0).getPosition());
            orientf[i].ME(pos.position(iMol));
            orientf[i].normalize();
            angleDisp += orientf[i].Mv1Squared(orient0[i]);
//            System.out.println("     "+Math.sqrt(orientf[i].Mv1Squared(orient0[i])));
        }
        
        disp = Math.sqrt(disp);
        angleDisp = Math.sqrt(angleDisp);
        System.out.println("\ndisp "+disp+"  angleDisp "+angleDisp);
           
        for(int i=0; i<sim.u.length; i++){
			sim.u[i] += deviation[i];
		}
        
        double newLatticeEnergy = sim.getEnergy(sim.u);
        System.out.println("\ndensity: " + density);
        System.out.println("Old Lattice Energy (per molecule): "+sim.initialLatticeEnergy);
        System.out.println("New Lattice Energy (per molecule): "+newLatticeEnergy);
        
        for(int i=0; i<sim.u.length; i++){
			if(i%5==0) System.out.print("{");
			System.out.print(sim.u[i]);
			if(i%5==4) {
				System.out.println("},");
			} else {
				System.out.print(", ");
			}
		}

        System.exit(1);
//        for (int l=0; l<201; l++) {
//            for (int i=0; i<4; i++) {
//                for (int k=i%2; k<numMolecule; k+=2) {
//                    boolean isA = (k/nA)%2 == 0;
//                    if ((i<2 && !isA) || (i>1 && isA)) continue;
//                    IMolecule iMol = sim.box.getMoleculeList().getMolecule(k);
//                    p.E(pos.position(iMol));
//                    for (int j=0; j<3; j++) {
//                        if (l==0) {
//                            p.setX(j, p.getX(j)-(xf[i*3+j] - x0[i*3+j])*10);
//                        }
//                        else {
//                            p.setX(j, p.getX(j)+(xf[i*3+j] - x0[i*3+j])/10.0);
//                        }
//                    }
//                    translator.setDestination(p);
//                    translator.actionPerformed(iMol);
//                }
//            }
//            double u = 0;
//            for (int i=0; i<4; i++) {
//                for (int k=(4*nA)+i%2; k<numMolecule; k+=2) {
//                    boolean isA = (k/nA)%2 == 0;
//                    if ((i<2 && !isA) || (i>1 && isA)) continue;
//                    IMolecule iMol = sim.box.getMoleculeList().getMolecule(k);
//                    meterPotentialEnergy.setTarget(iMol);
//                    u += meterPotentialEnergy.getDataAsScalar();
//                    break;
//                }
//            }
////            meterPotentialEnergy.setTarget((IMolecule)null);
////            System.out.println((l-100)/100.0+" "+(meterPotentialEnergy.getDataAsScalar()-newLatticeEnergy)/numMolecule);
//            System.out.println((l-100)/100.0+" "+(u/8-newLatticeEnergy/numMolecule));
//        }            
        System.exit(1);

	}
	
    protected static void doTransform(IMolecule molecule, IMoleculePositionDefinition posDef, Tensor rotationTensor) {
        IAtomList childList = molecule.getChildList();
        Vector com = posDef.position(molecule);
        for (int iChild = 0; iChild<childList.getAtomCount(); iChild++) {
            IAtom a = childList.getAtom(iChild);
            Vector r = a.getPosition();
            r.ME(com);
            rotationTensor.transform(r);
            r.PE(com);
        }
    }


	protected Box box, ghostBox;
	protected Space space;
	protected P2Nitrogen potential;
	protected CoordinateDefinitionNitrogen coordinateDef;
	protected Primitive primitive;
	protected SpeciesN2 species;
	protected Basis basis;
	protected double[] u;
	
	protected AtomActionTranslateBy translateBy;
	protected MoleculeChildAtomAction atomGroupActionTranslate;
	protected Vector lsPosition;
	protected double xVecBox, yVecBox, zVecBox, rC, density;  
	protected Vector[] boxVec;
	protected int nLayer;
    protected MoleculeChildAtomAction atomGroupAction;
    protected double initialLatticeEnergy;
	
	private static final long serialVersionUID = 1L;
}
