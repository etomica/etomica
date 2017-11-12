/**
 * Interfaces and classes that define elementary actions that can be performed
 * on simulation elements. Examples include classes that move an atom, or that change the density of a box.
 * This package defines the {@link etomica.action.IAction Action} interface. This interface includes the <tt>actionPerformed</tt>
 * method, which a class implements to define its action.
 * <p>
 * Interfaces extending Action define actions for different <i>etomica</i> elements, such as
 * Box, Atom, and so on.
 * <p>
 * {@link etomica.action.Activity Activity} implements Action and is an abstract class appropriate
 * to encapsulate more complex, time-consuming tasks.  Classes extending Activity are collected
 * in the subpackage {@link etomica.action.activity etomica.action.activity}.
 */
package etomica.action;
