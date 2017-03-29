/**
 * Provides interfaces and classes that define elementary actions that can be performed
 * on simulation elements. Examples include classes that move an atom, or that change the density of a box.
 * This package defines the {@link etomica.action.IAction Action} interface. This interface includes the <tt>actionPerformed</tt>
 * method, which a class implements to define its action.
 * <p>
 * Interfaces extending Action define actions for different <i>etomica</i> elements, such as
 * Box, Atom, and so on.  Atomset actions are expected as arguments
 * to the allAtoms methods of {@link etomica.atom.iterator atom iterators}.
 * <p>
 * {@link etomica.action.Activity Activity} implements Action and is an abstract class appropriate
 * to encapsulate more complex, time-consuming tasks.  Classes extending Activity are collected
 * in the subpackage {@link etomica.action.activity etomica.action.activity}.
 * <!--
 * <h2>Package Specification</h2>
 * ##### FILL IN ANY SPECS NEEDED BY JAVA COMPATIBILITY KIT #####
 * <ul>
 * <li><a href="">##### REFER TO ANY FRAMEMAKER SPECIFICATION HERE #####</a>
 * </ul>
 * <h2>Related Documentation</h2>
 * For overviews, tutorials, examples, guides, and tool documentation, please see:
 * <ul>
 * <li><a href="">##### REFER TO NON-SPEC DOCUMENTATION HERE #####</a>
 * </ul>
 * -->
 * <!-- Put @see and @since tags down here. -->
 */
package etomica.action;
