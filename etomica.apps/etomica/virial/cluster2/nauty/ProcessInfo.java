/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.nauty;

import java.util.List;

import etomica.virial.cluster2.graph.Tagged;

public interface ProcessInfo extends Tagged {

  /**
   * @return path where the process should run from.
   */
  public String getPath();

  /**
   * @return concatenation of the program and parameters for use by the
   *         ProcessWrapper.
   */
  public List<String> getCommand();

  /**
   * @return true if the error stream is being redirected to the output stream
   */
  public boolean redirectErrorStream();
}