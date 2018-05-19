/**
 * Define elementary actions that can be performed on simulation elements.
 * Examples include classes that move an Atom, or that change the density of a Box.
 * Classes defining such actions implement the {@link etomica.action.IAction Action} interface.
 * Interfaces extending IAction define actions for different simulation elements, such as
 * Box, Atom, and so on.
 * <p>
 * {@link etomica.action.Activity Activity} implements Action and is an abstract class appropriate
 * to encapsulate more complex, time-consuming tasks.  Classes extending Activity are collected
 * in the subpackage {@link etomica.action.activity etomica.action.activity}. Whereas Actions must
 * be completed before further processing continues, Activity instances generally operate on a different
 * thread than the one that calls them, and can be paused, resumed or halted before they complete.
 * Classes extending Activity are collected in the subpackage {@link etomica.action.activity etomica.action.activity}.
 */
package etomica.action;

