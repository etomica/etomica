/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.util.random.IRandom;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.space3d.RotationTensor3D;
import etomica.space3d.Tensor3D;
import etomica.units.Degree;
import etomica.util.random.RandomNumberGenerator;

public class TestRotationVector {

	public TestRotationVector(){
		space = Space.getInstance(3);
		axis = new Vector[3];
		for (int i=0; i<3; i++){
			axis[i] = space.makeVector();
		}
		
		axis[0].E(new double[]{1.0, 0.0, 0.0});
		axis[1].E(new double[]{0.0, 1.0, 0.0});
		axis[2].E(new double[]{0.0, 0.0, 1.0});
		
		u = new double[2];
		
		tensor = new Tensor3D(new double[][]{{1.0, 0.0, 0.0},{0.0, 1.0, 0.0},{0.0, 0.0, 1.0}});
		rotation = new RotationTensor3D();
		rotation.E(tensor);
		rotationAxis = space.makeVector();
	}
	
	public double[] calcU(Vector vector){
		System.out.println("In calcU");
		double u3 = vector.dot(axis[1]);
		double u4 = vector.dot(axis[2]);
		double ratio = Math.abs(u3/u4);
		
		double a = vector.dot(axis[0]);
		double theta = Math.acos(a);
		System.out.println("theta: " + Degree.UNIT.fromSim(theta));
		
		if(Degree.UNIT.fromSim(theta) > 179.999999){
			System.out.println("%%%%%%% theta is greater than 180 deg");
			u[0] = Math.sqrt(2);
			u[1] = Math.sqrt(2);
			for (int i=0; i<u.length; i++){
				System.out.println("CalcU u["+i+"]: "+u[i]);
			}
			System.out.println("END of calcU\n");
			return u;
		}
		   if(Math.abs(u4) > -1e-10 && Math.abs(u4) < 1e-10){
               u[0] = Math.sqrt(2*(1-Math.cos(theta)));
               if(u3 <0.0){
            	   u[0] = -u[0];
               }
               u[1] = u4;
		   } else {
               if(u4 < 0.0){
                       u[1] = -Math.sqrt(2*(1-Math.cos(theta))/(ratio*ratio+1));
               } else {
                       u[1] = Math.sqrt(2*(1-Math.cos(theta))/(ratio*ratio+1));
               }

               if (u3 < 0.0){
                       u[0] = -ratio*Math.sqrt(2*(1-Math.cos(theta))/(ratio*ratio+1));
               } else {
                       u[0] = ratio*Math.sqrt(2*(1-Math.cos(theta))/(ratio*ratio+1));
               }
       }

		for (int i=0; i<u.length; i++){
			System.out.println("CalcU u["+i+"]: "+u[i]);
		}
		System.out.println("END of calcU\n");
		return u; 
	}
	
	public void setToU(double[] u, Vector r){
		System.out.println("In setToU");
		double angle = Math.acos(r.dot(axis[0]));
		
		if(Math.sqrt(angle) > 1e-7){
			rotationAxis.E(r);
			rotationAxis.XE(axis[0]);
			rotationAxis.normalize();
			
			System.out.println("angle: " + angle);
			
			rotation.setRotationAxis(rotationAxis, angle);
			rotation.transform(r);
			r.normalize();
			System.out.println("Back to initial nominal position: "+ r.toString());
		}
		/*
		 * 
		 * 
		 * scale u3 and u4 accordingly so that they will satisfy the
		 *  condition u3^2 + u4^2 < 4.0
		 * 
		 */
		if(Math.abs(u[0]) > 1e-7 || Math.abs(u[1]) > 1e-7){
			double u0 = u[0];
			double u1 = u[1];
			double check = u0*u0+u1*u1;
			if((Math.abs(u0) > (Math.sqrt(2)+1e-10) || Math.abs(u1) > (Math.sqrt(2)+1e-10)) 
					&& (check > 3.99999999)){
				System.out.println("*****Free Rotor******");
				
				IRandom random = new RandomNumberGenerator();
				double randU0 = random.nextDouble()*Math.sqrt(2);
				double randU1 = random.nextDouble()*Math.sqrt(2);
				
				u0 = randU0;
				u1 = randU1;
				if (u[1] < 0.0){
					u[1] = -u1;
				} else {
					u[1] = u1;
				}
				if (u[0] < 0.0){
					u[0] = -u0;
				} else {
					u[0] = u0;
				}
				
			}
			
			System.out.println("u[]:" + u[0] + " " +u[1]);
			System.out.println("Math.acos component: " + (1 - (u[0]*u[0] + u[1]*u[1])*0.5));
			
			double theta = Math.acos(1.0000000000000004 - (u[0]*u[0] + u[1]*u[1])*0.5);
			System.out.println("theta in setToU: " + Degree.UNIT.fromSim(theta));
			
			r.E(0.0);
			r.PEa1Tv1(u[0], axis[1]);
			r.PEa1Tv1(u[1], axis[2]);
			r.normalize();
			
			rotationAxis.E(r);
			rotationAxis.XE(axis[0]);
			rotationAxis.normalize();
			rotation.setRotationAxis(rotationAxis, (Math.PI/2-theta));
			rotation.transform(r);		
		}
		System.out.println("END of setToU");
		System.out.println("After setToU destination: " + r.toString()+"\n");
	}
	
	public static void main (String[] args){
		TestRotationVector testVector = new TestRotationVector();
		
		Vector rVector = testVector.space.makeVector();
		rVector.E(new double[]{-0.99, 0.01, 0.01});
		rVector.normalize();
		System.out.println("Initial position: " + rVector.toString());
		double[] u = testVector.calcU(rVector);

		System.out.println("check: " + (u[0]*u[0]+u[1]*u[1]));
		for(int i=0; i<u.length; i++){
			u[i]/=1.00884;
		}
		//u = new double[]{-1.95, -0.44};
		System.out.println("sum u^2: " + (u[0]*u[0]+u[1]*u[1]));
		testVector.setToU(u, rVector);

		u = testVector.calcU(rVector);
		
	}
	
	protected double[] u;
	protected Vector[] axis;
	protected Space space;
	protected Tensor3D tensor;
	protected RotationTensor3D rotation;
	protected Vector rotationAxis;
	
	
	
}
