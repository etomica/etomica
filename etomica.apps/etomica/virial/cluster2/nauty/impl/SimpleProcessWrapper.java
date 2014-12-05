/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.nauty.impl;

import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.Map;

import etomica.virial.cluster2.nauty.ProcessInfo;
import etomica.virial.cluster2.nauty.ProcessWrapper;

public class SimpleProcessWrapper implements ProcessWrapper {

  private Reader output;
  private ProcessBuilder builder;
  private Map<String, String> environment;
  private ProcessInfo processInfo;

  public SimpleProcessWrapper(ProcessInfo info) {

    processInfo = info;
    builder = new ProcessBuilder(info.getCommand());
    environment = builder.environment();
  }

  public Map<String, String> getEnvironment() {

    return environment;
  }

  public Reader getProcessOutput() {

    return output;
  }

// @Override
// public List<String> getStdout() {
//
// return stdout;
// }
  public ProcessInfo getProcessInfo() {

    return processInfo;
  }

  /**
   * Blocking execution of the command.
   * 
   * @throws IOException
   */
  public void run() throws IOException {

    builder.redirectErrorStream(processInfo.redirectErrorStream());
    builder.directory(new File(processInfo.getPath()));
    Process process = builder.start();
    output = new InputStreamReader(process.getInputStream());
  }
}