/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.crystalviewer;

import etomica.space.Vector;
import junit.framework.TestCase;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.lattice.crystal.PrimitiveTetragonal;

public class BLCPrimitiveTetragonalLatticePlaneTest extends TestCase {

	private final int DEFAULT_SIZE = 5;
	private final int DEFAULT_MILLER[] = {0,0,1};
	private final int DEFAULT_BOX[] = {DEFAULT_SIZE, DEFAULT_SIZE, DEFAULT_SIZE};
	
	private String funcName = "";

	private double epsilon = 1.0E-5;;

	private LatticePlaneTestUtility lptu = null;

	public BLCPrimitiveTetragonalLatticePlaneTest(String name) {
		super(name);
		funcName = name;
	}

	protected void setUp() throws Exception {
		super.setUp();
		if (lptu == null) {
			lptu = new LatticePlaneTestUtility();			
	        lptu.createLatticeAndBox(lptu.TETRAGONAL, DEFAULT_MILLER, DEFAULT_BOX);
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
     * Miller indices = 0, 0, 1
     * size of cell (A, B, C) = 1.0
     * cells per side = 5
     * plane = 0.0
     */
    public void testStandard() {

    	int idx = 0;
    	double cubicSize = 1.0;
    	double plane = 0.0;
    	IAtomList leafList = null;

        lptu.createLatticeAndBox(lptu.TETRAGONAL, DEFAULT_MILLER, DEFAULT_BOX);

        ((PrimitiveTetragonal)lptu.getLattice().getPrimitive()).setSizeAB(cubicSize);
        ((PrimitiveTetragonal)lptu.getLattice().getPrimitive()).setSizeC(cubicSize);
        lptu.setDimensions(DEFAULT_SIZE);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());
        double spacePos = lptu.getLatticePlaneSpacePosition();

    	leafList = lptu.getBox().getLeafList();

    	try {
		    for(idx = 0; idx < leafList.size(); idx++) {
			    IAtom a =  leafList.get(idx);
                if(a.getPosition().getX(2) >= spacePos-epsilon &&
                   a.getPosition().getX(2) <= spacePos+epsilon) {
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
            if(a.getPosition().getX(2) >= spacePos-epsilon &&
               a.getPosition().getX(2) <= spacePos+epsilon) {
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
     * Miller indices = 0, 0, 1
     * size of cell (A, B) = 1.63
     * size of cell (C) = 1.0
     * cells per side = 5
     * plane = 0.0
     */
    public void testABCellSizeIncrease() {

    	int idx = 0;
    	double cubicSizeAB = 1.63;
    	double cubicSizeC = 1.0;
    	double plane = 0.0;
    	IAtomList leafList = null;

        lptu.createLatticeAndBox(lptu.TETRAGONAL, DEFAULT_MILLER, DEFAULT_BOX);

        ((PrimitiveTetragonal)lptu.getLattice().getPrimitive()).setSizeAB(cubicSizeAB);
        ((PrimitiveTetragonal)lptu.getLattice().getPrimitive()).setSizeC(cubicSizeC);
        lptu.setDimensions(DEFAULT_SIZE);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());
        double spacePos = lptu.getLatticePlaneSpacePosition();

    	leafList = lptu.getBox().getLeafList();

    	try {
		    for(idx = 0; idx < leafList.size(); idx++) {
			    IAtom a =  leafList.get(idx);
                if(a.getPosition().getX(2) >= spacePos-epsilon &&
                   a.getPosition().getX(2) <= spacePos+epsilon) {
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
            if(a.getPosition().getX(2) >= spacePos-epsilon &&
               a.getPosition().getX(2) <= spacePos+epsilon) {
            	System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should be in plane but is not.");
            }
            else {
            	System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should not be in plane but is.");       	
            }
         	fail();
        }
    } // End testABCellSizeIncrease()

    /*
     * Miller indices = 1, 3, 2
     * size of cell (A, B, C) = 1.0
     * cells per side = 9
     * plane = 7.0
     */
    public void testOddMillerIndicesDistantPlane() {

    	int idx = 0;
    	double cubicSize = 1.0;
    	double plane = 7.0;
    	IAtomList leafList = null;
    	int size = 9;
    	int itemsFound = 0;
    	int[] millerIndices = new int[] { 1, 3, 2 };
        double actualPlane[][] =
                  { { -4.0, 1.0, 4.0 }, { -4.0, 3.0, 1.0 }, { -3.0, 2.0, 2.0 }, 
                    { -3.0, 4.0, -1.0 }, { -2.0, 1.0, 3.0 }, { -2.0, 3.0, 0.0 }, 
                    { -1.0, 0.0, 4.0 }, { -1.0, 2.0, 1.0 }, { -1.0, 4.0, -2.0 }, 
                    { 0.0, 1.0, 2.0 }, { 0.0, 3.0, -1.0 }, { 1.0, 0.0, 3.0 }, 
                    { 1.0, 2.0, 0.0 }, { 1.0, 4.0, -3.0 }, { 2.0, -1.0, 4.0 }, 
                    { 2.0, 1.0, 1.0 }, { 2.0, 3.0, -2.0 }, { 3.0, 0.0, 2.0 }, 
                    { 3.0, 2.0, -1.0 }, { 3.0, 4.0, -4.0 }, { 4.0, -1.0, 3.0 }, 
                    { 4.0, 1.0, 0.0 }, { 4.0, 3.0, -3.0 } }; 

        DoubleTwoDArray dd = new DoubleTwoDArray(actualPlane);

        lptu.createLatticeAndBox(lptu.TETRAGONAL, millerIndices, new int[] {size, size, size});
        
        ((PrimitiveTetragonal)lptu.getLattice().getPrimitive()).setSizeAB(cubicSize);
        ((PrimitiveTetragonal)lptu.getLattice().getPrimitive()).setSizeC(cubicSize);
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

        assertEquals(actualPlane.length, itemsFound);

    } // End testCellSizeIncrease()

    /*
     * Miller indices = 0, 0, 1
     * size of cell (A, B, C) = 1.0
     * cells per side = 8
     * plane = 0.0
     */
    public void testEvenAtomsPerSideZeroPlane() {

    	int idx = 0;
    	double cubicSize = 1.0;
    	double plane = 0.0;
    	IAtomList leafList = null;
    	int dimensionSize = 8;

        lptu.createLatticeAndBox(lptu.TETRAGONAL, DEFAULT_MILLER, DEFAULT_BOX);

        ((PrimitiveTetragonal)lptu.getLattice().getPrimitive()).setSizeAB(cubicSize);
        ((PrimitiveTetragonal)lptu.getLattice().getPrimitive()).setSizeC(cubicSize);
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

        lptu.createLatticeAndBox(lptu.TETRAGONAL, millerIndices, DEFAULT_BOX);

        ((PrimitiveTetragonal)lptu.getLattice().getPrimitive()).setSizeAB(cubicSize);
        ((PrimitiveTetragonal)lptu.getLattice().getPrimitive()).setSizeC(cubicSize);
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
     * Miller indices = 1, 1, 2
     * size of cell (A, B) = 1.2
     * size of cell (C) = 1.4
     * cells per side = 6
     * plane = 1.95
     */
    public void testPlaneMinusFiveHundreths() {

    	int idx = 0;
    	double cubicSizeAB = 1.2;
    	double cubicSizeC = 1.4;
    	double plane = 1.95;
    	IAtomList leafList = null;
    	int size = 6;
    	int[] millerIndices = new int[] { 1, 1, 2 };

        lptu.createLatticeAndBox(lptu.TETRAGONAL, millerIndices, new int[] {size, size, size});
        
        ((PrimitiveTetragonal)lptu.getLattice().getPrimitive()).setSizeAB(cubicSizeAB);
        ((PrimitiveTetragonal)lptu.getLattice().getPrimitive()).setSizeC(cubicSizeC);
        lptu.setDimensions(size);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());

    	leafList = lptu.getBox().getLeafList();

    	try {
		    for(idx = 0; idx < leafList.size(); idx++) {
			    IAtom a = leafList.get(idx);

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

    } // End testPlaneMinusFiveHundreths()

    /*
     * Miller indices = 1, 1, 2
     * size of cell (A, B) = 1.2
     * size of cell (C) = 1.4
     * cells per side = 6
     * plane = 2.0
     */
    public void testPlane() {

    	int idx = 0;
    	double cubicSizeAB = 1.2;
    	double cubicSizeC = 1.4;
    	double plane = 2.0;
    	IAtomList leafList = null;
    	int size = 6;
    	int itemsFound = 0;
    	int[] millerIndices = new int[] { 1, 1, 2 };
        double actualPlane[][] =
                  { { -3.0, -0.6, 3.5 }, { -3.0, 1.8, 2.1}, { -1.8, -1.8, 3.5 },
        		    { -1.8, 0.6, 2.1}, { -1.8, 3.0, 0.7 }, { -0.6, -3.0, 3.5},
                    { -0.6, -0.6, 2.1 }, { -0.6, 1.8, 0.7}, { 0.6, -1.8, 2.1 },
                    { 0.6, 0.6, 0.7}, { 0.6, 3.0, -0.7 }, { 1.8, -3.0, 2.1},
                    { 1.8, -0.6, 0.7 }, { 1.8, 1.8, -0.7}, { 3.0, -1.8, 0.7 },
                    { 3.0, 0.6, -0.7}, { 3.0, 3.0, -2.1} };

        DoubleTwoDArray dd = new DoubleTwoDArray(actualPlane);

        lptu.createLatticeAndBox(lptu.TETRAGONAL, millerIndices, new int[] {size, size, size});
        
        ((PrimitiveTetragonal)lptu.getLattice().getPrimitive()).setSizeAB(cubicSizeAB);
        ((PrimitiveTetragonal)lptu.getLattice().getPrimitive()).setSizeC(cubicSizeC);
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

        assertEquals(actualPlane.length, itemsFound);

    } // End testPlane()

    /*
     * Miller indices = 1, 1, 2
     * size of cell (A, B) = 1.2
     * size of cell (C) = 1.4
     * cells per side = 6
     * plane = 2.05
     */
    public void testPlanePlusFiveHundreths() {

    	int idx = 0;
    	double cubicSizeAB = 1.2;
    	double cubicSizeC = 1.4;
    	double plane = 2.05;
    	IAtomList leafList = null;
    	int size = 6;
    	int[] millerIndices = new int[] { 1, 1, 2 };

        lptu.createLatticeAndBox(lptu.TETRAGONAL, millerIndices, new int[] {size, size, size});
        
        ((PrimitiveTetragonal)lptu.getLattice().getPrimitive()).setSizeAB(cubicSizeAB);
        ((PrimitiveTetragonal)lptu.getLattice().getPrimitive()).setSizeC(cubicSizeC);
        lptu.setDimensions(size);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());

    	leafList = lptu.getBox().getLeafList();

    	try {
		    for(idx = 0; idx < leafList.size(); idx++) {
			    IAtom a = leafList.get(idx);

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

    } // End testPlanePlusFiveHundreths()

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
