/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.crystalviewer;

import etomica.space.Vector;
import junit.framework.TestCase;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.lattice.BravaisLattice;
import etomica.lattice.crystal.PrimitiveCubic;


public class LatticeCubicSimpleLatticePlaneTest extends TestCase {

	BravaisLattice[] lattices = null;


	private final int DEFAULT_SIZE = 7;
	private final int DEFAULT_MILLER[] = {1,0,0};
	private final int DEFAULT_BOX[] = {DEFAULT_SIZE, DEFAULT_SIZE, DEFAULT_SIZE};
	
	private String funcName = "";

	private double epsilon = 1.0E-5;;

	private LatticePlaneTestUtility lptu = null;

	public LatticeCubicSimpleLatticePlaneTest(String name) {
		super(name);
		funcName = name;
	}

	protected void setUp() throws Exception {
    	
    	super.setUp();

		if (lptu == null) {
			lptu = new LatticePlaneTestUtility();			
	        lptu.createLatticeAndBox(lptu.SIMPLE_CUBIC, DEFAULT_MILLER, DEFAULT_BOX);
	        lptu.setDimensions(DEFAULT_SIZE);
		}

	}

	protected void tearDown() throws Exception {
		super.tearDown();
	}

	private double[] makeArray(Vector v) {
	    return new double[] {v.getX(0), v.getX(1), v.getX(2)};
	}

    /*
     * Miller indices = 1, 0, 0
     * size of cell (A, B, C) = 1.0
     * cells per side = 7
     * plane = 2.0
     */
    public void testStandard() {

    	int idx = 0;
    	double cubicSize = 1.0;
    	double plane = 2.0;
    	IAtomList leafList = null;

        lptu.createLatticeAndBox(lptu.SIMPLE_CUBIC, DEFAULT_MILLER, DEFAULT_BOX);

        ((PrimitiveCubic)lptu.getLattice().getPrimitive()).setSizeABC(cubicSize);

        lptu.setDimensions(DEFAULT_SIZE);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());
        double spacePos = lptu.getLatticePlaneSpacePosition();

        leafList = lptu.getBox().getLeafList();

    	try {
		    for(idx = 0; idx < leafList.size(); idx++) {
			    IAtom a =  leafList.get(idx);
                if(a.getPosition().getX(0) >= spacePos-epsilon &&
                   a.getPosition().getX(0) <= spacePos+epsilon) {
            	    assertTrue(lptu.getLatticePlane().inPlane(
            	    		a.getPosition()));
                }
                else {
            	    assertFalse(lptu.getLatticePlane().inPlane(
            	    		a.getPosition()));
                }
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtom a =  leafList.get(idx);
            if(a.getPosition().getX(0) >= spacePos-epsilon &&
               a.getPosition().getX(0) <= spacePos+epsilon) {
            	System.out.println(funcName + " -> Atom position : " + a.getPosition() +
            			" should be in plane but is not.");
            }
            else {
            	System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should not be in plane but is.");       	
            }
         	fail();
        }

    } // End testStandard()

    /*
     * Miller indices = 1, 0, 0
     * size of cell (A, B, C) = 1.63
     * cells per side = 7
     * plane = 1.0
     */
    public void testCellSizeIncrease() {

    	int idx = 0;
    	double cubicSize = 1.63;
    	double plane = 1.0;
    	IAtomList leafList = null;

        lptu.createLatticeAndBox(lptu.SIMPLE_CUBIC, DEFAULT_MILLER, DEFAULT_BOX);

        ((PrimitiveCubic)lptu.getLattice().getPrimitive()).setSizeABC(cubicSize);

        lptu.setDimensions(DEFAULT_SIZE);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());
        double spacePos = lptu.getLatticePlaneSpacePosition();

        leafList = lptu.getBox().getLeafList();

    	try {
		    for(idx = 0; idx < leafList.size(); idx++) {
			    IAtom a =  leafList.get(idx);
                if(a.getPosition().getX(0) >= spacePos-epsilon &&
                   a.getPosition().getX(0) <= spacePos+epsilon) {
            	    assertTrue(lptu.getLatticePlane().inPlane(
            	    		a.getPosition()));
                }
                else {
            	    assertFalse(lptu.getLatticePlane().inPlane(
            	    		a.getPosition()));
                }
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtom a =  leafList.get(idx);
            if(a.getPosition().getX(0) >= spacePos-epsilon &&
               a.getPosition().getX(0) <= spacePos+epsilon) {
            	System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should be in plane but is not.");
            }
            else {
            	System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should not be in plane but is.");       	
            }
         	fail();
        }
    } // End testCellSizeIncrease()

    /*
     * Miller indices = 1, 3, 2
     * size of cell (A, B, C) = 1.0
     * cells per side = 10
     * plane = 7.0
     */
    public void testOddMillerIndicesDistantPlane() {

    	int idx = 0;
    	double cubicSize = 1.0;
    	double plane = 7.0;
    	IAtomList leafList = null;
    	int size = 10;
    	int itemsFound = 0;
    	int[] millerIndices = new int[] { 1, 3, 2 };
        double actualPlane[][] =
              { {-4.5, 1.5, 3.5}, {-4.5, 3.5, 0.5}, {-3.5, 0.5, 4.5},
        		{ -3.5, 2.5, 1.5},{ -3.5, 4.5, -1.5},{ -2.5, 1.5, 2.5},
        		{ -2.5, 3.5, -0.5},{ -1.5, 0.5, 3.5},{ -1.5, 2.5, 0.5},
        		{ -1.5, 4.5, -2.5},{ -0.5, -0.5, 4.5},{ -0.5, 1.5, 1.5},
        		{ -0.5, 3.5, -1.5},{ 0.5, 0.5, 2.5},{ 0.5, 2.5, -0.5},
        		{ 0.5, 4.5, -3.5},{ 1.5, -0.5, 3.5},{ 1.5, 1.5, 0.5},
        		{ 1.5, 3.5, -2.5},{ 2.5, -1.5, 4.5},{ 2.5, 0.5, 1.5},
        		{ 2.5, 2.5, -1.5},{ 2.5, 4.5, -4.5},{ 3.5, -0.5, 2.5},
        		{ 3.5, 1.5, -0.5},{ 3.5, 3.5, -3.5},{ 4.5, -1.5, 3.5},
        		{ 4.5, 0.5, 0.5},{ 4.5, 2.5, -2.5} };

        DoubleTwoDArray dd = new DoubleTwoDArray(actualPlane);

        lptu.createLatticeAndBox(lptu.SIMPLE_CUBIC, millerIndices, new int[] {size, size, size});
        
        ((PrimitiveCubic)lptu.getLattice().getPrimitive()).setSizeABC(cubicSize);

        lptu.setDimensions(size);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());

    	leafList = lptu.getBox().getLeafList();

    	try {
		    for(idx = 0; idx < leafList.size(); idx++) {
			    IAtom a = leafList.get(idx);

			    if(dd.contains(makeArray(a.getPosition())) == true) {
			    	itemsFound++;
            	    assertTrue(lptu.getLatticePlane().inPlane(
            	    		a.getPosition()));
                }
                else {
            	    assertFalse(lptu.getLatticePlane().inPlane(
            	    		a.getPosition()));
                }
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtom a =  leafList.get(idx);
            if(dd.contains(makeArray(a.getPosition()))) {
            	System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should be in plane but is not.");
            }
            else {
            	System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should not be in plane but is.");       	
            }
         	fail();
        }

        assertEquals(itemsFound, actualPlane.length);

    } // End testCellSizeIncrease()

    /*
     * Miller indices = 1, 0, 0
     * size of cell (A, B, C) = 1.0
     * cells per side = 8
     * plane = 1.0
     */
    public void testEvenAtomsPerSideZeroPlane() {

    	int idx = 0;
    	double cubicSize = 1.0;
    	double plane = 1.0;
    	IAtomList leafList = null;
    	int dimensionSize = 8;

        lptu.createLatticeAndBox(lptu.SIMPLE_CUBIC, DEFAULT_MILLER,
        		                   new int[] {dimensionSize, dimensionSize, dimensionSize});

        ((PrimitiveCubic)lptu.getLattice().getPrimitive()).setSizeABC(cubicSize);

        lptu.setDimensions(dimensionSize);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());

    	leafList = lptu.getBox().getLeafList();

    	try {
		    for(idx = 0; idx < leafList.size(); idx++) {
			    IAtom a =  leafList.get(idx);
            	assertFalse(lptu.getLatticePlane().inPlane(
            	    	a.getPosition()));
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtom a =  leafList.get(idx);
            System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should not be in plane but is.");
         	fail();
        }
    	
    } // End testEvenAtomsPerSideZeroPlane

    /*
     * Miller indices = 0, 1, 0
     * size of cell (A, B, C) = 1.0
     * cells per side = 8
     * plane = 0.5
     */
    public void testEvenAtomsPerSideZeroPt5Plane() {

    	int idx = 0;
    	double cubicSize = 1.0;
    	double plane = 0.5;
    	IAtomList leafList = null;
    	int dimensionSize = 8;
    	int[] millerIndices = new int[] { 0, 1, 0 };

        lptu.createLatticeAndBox(lptu.SIMPLE_CUBIC, millerIndices,
        		                   new int[] {dimensionSize, dimensionSize, dimensionSize});

        ((PrimitiveCubic)lptu.getLattice().getPrimitive()).setSizeABC(cubicSize);
        lptu.setDimensions(dimensionSize);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());
        double spacePos = lptu.getLatticePlaneSpacePosition();

    	leafList = lptu.getBox().getLeafList();

    	try {
		    for(idx = 0; idx < leafList.size(); idx++) {
			    IAtom a = leafList.get(idx);

                if(a.getPosition().getX(1) >= spacePos-epsilon &&
                   a.getPosition().getX(1) <= spacePos+epsilon) {
            	    assertTrue(lptu.getLatticePlane().inPlane(
            	    		a.getPosition()));
                }
                else {
            	    assertFalse(lptu.getLatticePlane().inPlane(
            	    		a.getPosition()));
                }
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtom a =  leafList.get(idx);
            if(a.getPosition().getX(1) >= spacePos-epsilon &&
               a.getPosition().getX(1) <= spacePos+epsilon) {
                System.out.println(funcName + " -> Atom position : " + a.getPosition() +
                 			" should be in plane but is not.");
             }
             else {
                 System.out.println(funcName + " ->Atom position : " + a.getPosition() +
                 			" should not be in plane but is.");       	
             }
         	fail();
        }
    	
    } // End testEvenAtomsPerSideZeroPt5Plane

    /*
     * Miller indices = 4, 1, 3
     * size of cell (A, B, C) = 1.45
     * cells per side = 9
     * plane = 3.95
     */
    public void testPlaneMinusFiveHundreths() {

    	int idx = 0;
    	double cubicSize = 1.45;
    	double plane = 3.95;
    	int[] miller = { 4, 1, 3 };
    	IAtomList leafList = null;
    	int dimensionSize = 9;

        lptu.createLatticeAndBox(lptu.SIMPLE_CUBIC, miller,
        		                   new int[] {dimensionSize, dimensionSize, dimensionSize});

        ((PrimitiveCubic)lptu.getLattice().getPrimitive()).setSizeABC(cubicSize);
        lptu.setDimensions(dimensionSize);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());

    	leafList = lptu.getBox().getLeafList();

    	try {
		    for(idx = 0; idx < leafList.size(); idx++) {
			    IAtom a =  leafList.get(idx);
            	assertFalse(lptu.getLatticePlane().inPlane(
            	    	a.getPosition()));
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtom a =  leafList.get(idx);
            System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should not be in plane but is.");  
            fail();
        }
    	
    } // End testPlaneMinusOneHundreth

    /*
     * Miller indices = 4, 1, 3
     * size of cell (A, B, C) = 1.45
     * cells per side = 9
     * plane = 4.0
     */
    public void testPlane() {

    	int idx = 0;
    	double cubicSize = 1.45;
    	double plane = 4.0;
    	int[] miller = { 4, 1, 3 };
    	IAtomList leafList = null;
    	int dimensionSize = 9;
    	int itemsFound = 0;
        double actualPlane[][] =
           { { -4.35, 5.8, 5.8 }, { -2.9, 0.0, 5.8 },
             { -2.9, 4.35, 4.35 }, { -1.45, -5.8, 5.8 },
             { -1.45, -1.45, 4.35 }, { -1.45, 2.9, 2.9 },
             { 0.0, -2.9, 2.9 }, { 0.0, 1.45, 1.45 },
             { 0.0, 5.8, 0.0 }, { 1.45, -4.35, 1.45 },
             { 1.45, 0.0, 0.0 }, { 1.45, 4.35, -1.45 },
             { 2.9, -5.8, 0.0 }, { 2.9, -1.45, -1.45 },
             { 2.9, 2.9, -2.9 }, { 4.350, -2.9, -2.9 },
             { 4.35, 1.45, -4.35 }, { 4.35, 5.8, -5.8 },
             { 5.8, -4.35, -4.35 }, { 5.8, 0.0, -5.8 } };

        DoubleTwoDArray dd = new DoubleTwoDArray(actualPlane);

        lptu.createLatticeAndBox(lptu.SIMPLE_CUBIC, miller,
        		                   new int[] {dimensionSize, dimensionSize, dimensionSize});

        ((PrimitiveCubic)lptu.getLattice().getPrimitive()).setSizeABC(cubicSize);
        lptu.setDimensions(dimensionSize);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());

    	leafList = lptu.getBox().getLeafList();

    	try {
		    for(idx = 0; idx < leafList.size(); idx++) {
			    IAtom a = leafList.get(idx);

			    if(dd.contains(makeArray(a.getPosition())) == true) {
                    itemsFound++;
            	    assertTrue(lptu.getLatticePlane().inPlane(
            	    		a.getPosition()));
                }
                else {
                	assertFalse(lptu.getLatticePlane().inPlane(
            	    		a.getPosition()));
                }
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtom a =  leafList.get(idx);
            if(dd.contains(makeArray(a.getPosition()))) {
            	System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should be in plane but is not.");
            }
            else {
            	System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should not be in plane but is.");       	
            }
         	fail();
        }
        
        assertEquals(actualPlane.length, itemsFound);

    } // End testPlane

    /*
     * Miller indices = 4, 1, 3
     * size of cell (A, B, C) = 1.45
     * cells per side = 9
     * plane = 4.05
     */
    public void testPlanePlusFiveHundreths() {

    	int idx = 0;
    	double cubicSize = 1.45;
    	double plane = 4.05;
    	int[] miller = { 4, 1, 3 };
    	IAtomList leafList = null;
    	int dimensionSize = 9;

        lptu.createLatticeAndBox(lptu.SIMPLE_CUBIC, miller,
        		                   new int[] {dimensionSize, dimensionSize, dimensionSize});

        ((PrimitiveCubic)lptu.getLattice().getPrimitive()).setSizeABC(cubicSize);
        lptu.setDimensions(dimensionSize);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());

    	leafList = lptu.getBox().getLeafList();

    	try {
		    for(idx = 0; idx < leafList.size(); idx++) {
			    IAtom a =  leafList.get(idx);
            	assertFalse(lptu.getLatticePlane().inPlane(
            	    	a.getPosition()));
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtom a =  leafList.get(idx);
            System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should not be in plane but is.");
         	fail();
        }
    	
    } // End testPlanePlusOneHundreth

    public class DoubleTwoDArray {
    	private double[][] array;
    	private double epsilon = 1.0E-5;

    	public DoubleTwoDArray(double[][] in) {
    		array = in;
    	}
    	
    	public boolean contains(double[] val) {
    		boolean b = false;

    		if(array[0].length == val.length) {
    			for(int i = 0; i < array.length; i++) {
    				boolean matching = true;
    				for(int j = 0; j < array[0].length; j++) {
    					if(val[j] < array[i][j] - epsilon ||
    					   val[j] > array[i][j] + epsilon) {
    						matching = false;
    						break;
    					}
    				}
    				if(matching == true) {
    					b = true;
    					break;
    				}
    			}
    		}
            return b;
    	}
    	
    	public void setEpsilon(double e) {
    		this.epsilon = e;
    	}

    }  // End public class DoubleTwoDArray

}
