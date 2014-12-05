/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.mvc.launchers;


public class ApplicationLauncher {

  /**
   * A simple entry point to run the application. In the future, one might want
   * to use the arguments to the main method to modify the way the UI is
   * constructed.
   */
  public static void main(String args[]) {

    ViewFactory.showView(ViewFactory.VN_VIRIAL);
  }
}