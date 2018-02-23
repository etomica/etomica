/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.crystalviewer;

import etomica.space.Vector;
import junit.framework.TestCase;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.lattice.crystal.PrimitiveHexagonal;

public class BLCPrimitiveHexagonalLatticePlaneTest extends TestCase {

	private final int DEFAULT_SIZE = 5;
	private final int DEFAULT_MILLER[] = {0,1,0};
	private final int DEFAULT_BOX[] = {DEFAULT_SIZE, DEFAULT_SIZE, DEFAULT_SIZE};
	
	private String funcName = "";

	private double epsilon = 1.0E-5;

	private LatticePlaneTestUtility lptu = null;

	public BLCPrimitiveHexagonalLatticePlaneTest(String name) {
		super(name);
		funcName = name;
	}

	protected void setUp() throws Exception {
		super.setUp();
		if (lptu == null) {
			lptu = new LatticePlaneTestUtility();			
	        lptu.createLatticeAndBox(lptu.HEXAGONAL, DEFAULT_MILLER, DEFAULT_BOX);
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
     * Miller indices = 0, 1, 0
     * size of cell (A, B, C) = 1.0
     * cells per side = 5
     * plane = 1.0
     */
    public void testStandard() {

    	int idx = 0;
    	double cubicSize = 1.0;
    	double plane = 1.0;
    	IAtomList leafList = null;

        lptu.createLatticeAndBox(lptu.HEXAGONAL, DEFAULT_MILLER, DEFAULT_BOX);


        ((PrimitiveHexagonal)lptu.getLattice().getPrimitive()).setSizeAB(cubicSize);
        ((PrimitiveHexagonal)lptu.getLattice().getPrimitive()).setSizeC(cubicSize);
        lptu.setDimensions(DEFAULT_SIZE);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());
        double spacePos = lptu.getLatticePlaneSpacePosition();

        //Only (and all) atoms where x = 0 should be in plane.
    	leafList = lptu.getBox().getLeafList();

    	try {
		    for(idx = 0; idx < leafList.size(); idx++) {
			    IAtom a = leafList.get(idx);
                if(a.getPosition().getX(1) >= spacePos-epsilon &&
                   a.getPosition().getX(1) <= spacePos+epsilon) {
            	    assertTrue(lptu.getLatticePlane().inPlane(a.getPosition()));
                }
                else {
            	    assertFalse(lptu.getLatticePlane().inPlane(a.getPosition()));
                }
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtom a = leafList.get(idx);
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

    } // End testStandard()

    /*
     * Miller indices = 0, 1, 0
     * size of cell (A, B) = 1.0
     * size of cell (A, B) = 2.3
     * cells per side = 5
     * plane = -1.0
     */
    public void testCCellSizeIncrease() {

    	int idx = 0;
    	double cubicSizeAB = 1.0;
    	double cubicSizeC = 2.3;
    	double plane = -1.0;
    	IAtomList leafList = null;

        lptu.createLatticeAndBox(lptu.HEXAGONAL, DEFAULT_MILLER, DEFAULT_BOX);

        ((PrimitiveHexagonal)lptu.getLattice().getPrimitive()).setSizeAB(cubicSizeAB);
        ((PrimitiveHexagonal)lptu.getLattice().getPrimitive()).setSizeC(cubicSizeC);
        lptu.setDimensions(DEFAULT_SIZE);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());
        double spacePos = lptu.getLatticePlaneSpacePosition();

        leafList = lptu.getBox().getLeafList();

    	try {
		    for(idx = 0; idx < leafList.size(); idx++) {
			    IAtom a =  leafList.get(idx);
                if(a.getPosition().getX(1) >= spacePos-epsilon &&
                   a.getPosition().getX(1) <= spacePos+epsilon) {
            	    assertTrue(lptu.getLatticePlane().inPlane(a.getPosition()));
                }
                else {
            	    assertFalse(lptu.getLatticePlane().inPlane(a.getPosition()));
                }
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtom a =  leafList.get(idx);
            if(a.getPosition().getX(1) >= spacePos-epsilon &&
               a.getPosition().getX(1) <= spacePos+epsilon) {
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
    	IAtomList leafList = null;
    	int size = 9;
    	int itemsFound = 0;
    	double plane = 7.0;
    	int[] millerIndices = new int[] { 1, 3, 2 };
        double actualPlane[][] =
         { { -4.5, 0.8660254037844388, 4.0 },{ -5.5, 2.598076211353317, 1.0 },
           { -4.0, 1.7320508075688772, 2.0 },{ -5.0, 3.4641016151377557, -1.0 },
           { -2.5, 0.866025403784439, 3.0 },{ -3.5, 2.5980762113533173, 0.0 },
           { -1.0, 0.0, 4.0 },{ -2.0, 1.7320508075688772, 1.0 },
           { -3.0, 3.4641016151377557, -2.0 },{ -0.5, 0.8660254037844388, 2.0 },
           { -1.5, 2.598076211353317, -1.0 },{ 1.0, 0.0, 3.0 },
           { 0.0, 1.7320508075688772, 0.0 },{ -1.0, 3.4641016151377557, -3.0 },
           { 2.5, -0.8660254037844393, 4.0 },{ 1.5, 0.8660254037844388, 1.0 },
           { 0.5, 2.5980762113533173, -2.0 },{ 3.0, 0.0, 2.0 },
           { 2.0, 1.7320508075688772, -1.0 },{ 1.0, 3.4641016151377557, -4.0 },
           { 4.5, -0.8660254037844393, 3.0 },{ 3.5, 0.8660254037844388, 0.0 },
           { 2.5, 2.5980762113533173, -3.0 } };

        DoubleTwoDArray dd = new DoubleTwoDArray(actualPlane);

        lptu.createLatticeAndBox(lptu.HEXAGONAL, millerIndices, new int[] {size, size, size});

        ((PrimitiveHexagonal)lptu.getLattice().getPrimitive()).setSizeAB(cubicSize);
        ((PrimitiveHexagonal)lptu.getLattice().getPrimitive()).setSizeC(cubicSize);
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
            	    assertTrue(lptu.getLatticePlane().inPlane(a.getPosition()));
                }
                else {
            	    assertFalse(lptu.getLatticePlane().inPlane(a.getPosition()));
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
     * Miller indices = 0, 1, 0
     * size of cell (A, B, C) = 1.0
     * cells per side = 4
     * plane = 0.0
     */
    public void testEvenAtomsPerSideZeroPlane() {

    	int idx = 0;
    	double cubicSize = 1.0;
    	IAtomList leafList = null;
    	int dimensionSize = 4;
    	double plane = 0.0;

        lptu.createLatticeAndBox(lptu.HEXAGONAL, DEFAULT_MILLER,
        		                   new int[] {dimensionSize, dimensionSize, dimensionSize});

        ((PrimitiveHexagonal)lptu.getLattice().getPrimitive()).setSizeAB(cubicSize);
        ((PrimitiveHexagonal)lptu.getLattice().getPrimitive()).setSizeC(cubicSize);
        lptu.setDimensions(dimensionSize);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());

    	leafList = lptu.getBox().getLeafList();

    	try {
		    for(idx = 0; idx < leafList.size(); idx++) {
			    IAtom a =  leafList.get(idx);
            	assertFalse(lptu.getLatticePlane().inPlane(a.getPosition()));
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
     * cells per side = 4
     * plane = 0.5
     */
    public void testEvenAtomsPerSideZeroPt5Plane() {

    	int idx = 0;
    	double cubicSize = 1.0;
    	IAtomList leafList = null;
    	int dimensionSize = 4;
    	double plane = 0.5;

        lptu.createLatticeAndBox(lptu.HEXAGONAL, DEFAULT_MILLER,
        		                   new int[] {dimensionSize, dimensionSize, dimensionSize});

        ((PrimitiveHexagonal)lptu.getLattice().getPrimitive()).setSizeAB(cubicSize);
        ((PrimitiveHexagonal)lptu.getLattice().getPrimitive()).setSizeC(cubicSize);
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
            	    assertTrue(lptu.getLatticePlane().inPlane(a.getPosition()));
                }
                else {
            	    assertFalse(lptu.getLatticePlane().inPlane(a.getPosition()));
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
     * size of cell (A, B) = 1.4
     * size of cell (C) = 1.4
     * cells per side = 5
     * plane = 0.05
     */
    public void testPlaneMinusFiveHundreth() {

    	int idx = 0;
    	double cubicSizeAB = 1.4;
    	double cubicSizeC = 1.4;
    	IAtomList leafList = null;
    	int size = 5;
    	double plane = 0.05;
    	int[] millerIndices = new int[] { 1, 1, 2 };

        lptu.createLatticeAndBox(lptu.HEXAGONAL, millerIndices, new int[] {size, size, size});

        ((PrimitiveHexagonal)lptu.getLattice().getPrimitive()).setSizeAB(cubicSizeAB);
        ((PrimitiveHexagonal)lptu.getLattice().getPrimitive()).setSizeC(cubicSizeC);
        lptu.setDimensions(size);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());

        leafList = lptu.getBox().getLeafList();

    	try {
		    for(idx = 0; idx < leafList.size(); idx++) {
			    IAtom a = leafList.get(idx);

            	assertFalse(lptu.getLatticePlane().inPlane(a.getPosition()));
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtom a =  leafList.get(idx);
            System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should not be in plane but is.");
        	fail();
        }

    } // End testPlaneMinusFiveHundreth()

    /*
     * Miller indices = 1, 1, 2
     * size of cell (A, B) = 1.4
     * size of cell (C) = 1.4
     * cells per side = 5
     * plane = 0.0
     */
    public void testPlane() {

    	int idx = 0;
    	double cubicSizeAB = 1.4;
    	double cubicSizeC = 1.4;
    	IAtomList leafList = null;
    	int size = 5;
    	int itemsFound = 0;
    	double plane = 0.0;
    	int[] millerIndices = new int[] { 1, 1, 2 };
        double actualPlane[][] =
                { { -1.4, -2.4248711305964283, 2.8 }, { -2.8, 0.0, 1.4 },
                  { -4.2, 2.4248711305964283, 0.0 }, { -0.7, -1.2124355652982142, 1.4 },
                  { -2.1, 1.2124355652982142, 0.0 }, { 1.4, -2.4248711305964283, 1.4 },
                  { 0.0, 0.0, 0.0 }, { -1.4, 2.4248711305964283, -1.4 },
                  { 2.1, -1.2124355652982142, 0.0 }, { 0.7, 1.2124355652982142, -1.4 },
                  { 4.2, -2.4248711305964283, 0.0 }, { 2.8, 0.0, -1.4 },
                  { 1.4, 2.4248711305964283, -2.8 } };
        DoubleTwoDArray dd = new DoubleTwoDArray(actualPlane);

        lptu.createLatticeAndBox(lptu.HEXAGONAL, millerIndices, new int[] {size, size, size});

        ((PrimitiveHexagonal)lptu.getLattice().getPrimitive()).setSizeAB(cubicSizeAB);
        ((PrimitiveHexagonal)lptu.getLattice().getPrimitive()).setSizeC(cubicSizeC);
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
            	    assertTrue(lptu.getLatticePlane().inPlane(a.getPosition()));
                }
                else {
            	    assertFalse(lptu.getLatticePlane().inPlane(a.getPosition()));
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
     * size of cell (A, B) = 1.4
     * size of cell (C) = 1.4
     * cells per side = 5
     * plane = 0.05
     */
    public void testPlanePlusFiveHundreth() {

    	int idx = 0;
    	double cubicSizeAB = 1.4;
    	double cubicSizeC = 1.4;
    	IAtomList leafList = null;
    	int size = 5;
    	double plane = 0.05;
    	int[] millerIndices = new int[] { 1, 1, 2 };

        lptu.createLatticeAndBox(lptu.HEXAGONAL, millerIndices, new int[] {size, size, size});

        ((PrimitiveHexagonal)lptu.getLattice().getPrimitive()).setSizeAB(cubicSizeAB);
        ((PrimitiveHexagonal)lptu.getLattice().getPrimitive()).setSizeC(cubicSizeC);
        lptu.setDimensions(size);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());

        leafList = lptu.getBox().getLeafList();

    	try {
		    for(idx = 0; idx < leafList.size(); idx++) {
			    IAtom a = leafList.get(idx);

            	assertFalse(lptu.getLatticePlane().inPlane(a.getPosition()));
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtom a =  leafList.get(idx);
            System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should not be in plane but is.");
        	fail();
        }

    } // End testPlanePlusFiveHundreth()

    public static class DoubleTwoDArray {
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
