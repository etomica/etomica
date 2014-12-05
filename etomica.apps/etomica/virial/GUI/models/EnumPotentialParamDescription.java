/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI.models;

import etomica.units.Kelvin;

public enum EnumPotentialParamDescription {
	
	
	//Array Order so far
	/*
	 * 
	 * 
	 * LJ
	 * CO2
	 * C in EPM2
	 * O in EPM2
	 * CO in EPM2
	 * C in Trappe
	 * O in Trappe
	 * CO in Trappe
	 * 
	 * 
	 */
	SIGMA(0.0,"Potential Well-Depth", "("+Character.toString((char) 193)+")"),
	
	EPSILON(0.0, "Distance at which Interatomic potential between particles is zero","(K)"),
	
	MOMENT(0.0, "Moment","("+Character.toString((char) 283)+")"),
	
	BONDL(0.0, "Bond Length","("+Character.toString((char) 193)+")"),
	
	MOMENTSQR(0.0,"Moment square",null),
	
	CHARGE(0.0,"Charge","(e)"),
	
	NominalbondL(0.0,"Fixed bond length between neighboring pseudo-atoms","("+Character.toString((char) 193)+")"),
	
	theta(0.0,"Equilibrium bond angle","degrees"),
	
	forceconstant(0.0,"Force Constant(k0/kB)","K/rad2"),
	
	TEMPERATURE(0.0,"Temperature of Simulation","(K)"),
	
	STEPS(0.0,"No of Steps",null),
	
	SIGMAHSREF(0.0,"Hard Sphere Reference",null),
	
	NUMBER(0.0,"No of Spheres",null);
	

	private Double defaultValues;
	private String description;
	private String unit;
	
	EnumPotentialParamDescription(Double Value, String Description, String Unit){
		this.defaultValues = Value;
		this.description = Description;
		this.unit = Unit;
	}
	
	public Double defaultValue(){
		return defaultValues;
	}
	
	public String description(){
		return description;
	}
	
	public String unit(){
		return unit;
	}
}
