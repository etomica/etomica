/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.nauty;

import java.io.IOException;
import java.io.Reader;
import java.util.Map;

public interface ProcessWrapper {

  public Map<String, String> getEnvironment();

  public ProcessInfo getProcessInfo();

  public Reader getProcessOutput();

  public void run() throws IOException;
}