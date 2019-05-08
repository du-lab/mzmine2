/*
 * Copyright (C) 2018 Du-Lab Team <dulab.binf@gmail.com>
 *
 * This program is free software; you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
 * even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if
 * not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.
 */

package net.sf.mzmine.modules.peaklistmethods.dataanalysis.mummichog;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import driver.ExecuteMummiChog;
import net.sf.mzmine.datamodel.*;
import net.sf.mzmine.datamodel.impl.SimplePeakIdentity;
import net.sf.mzmine.taskcontrol.AbstractTask;
import net.sf.mzmine.taskcontrol.TaskStatus;
import pojo.Compound;
import pojo.MummichogParams;

public class MummichogTask extends AbstractTask {

  private static final String P_VALUE_KEY = "STUDENT_P_VALUE";
  private static final String T_VALUE_KEY = "STUDENT_T_VALUE";

  private Logger logger = Logger.getLogger(this.getClass().getName());
  private double finishedPercentage = 0.0;

  private final PeakListRow[] peakListRows;
  private final double cutoff;
  private final String network;
  private final boolean force_primary_ion;
  private final String modeling;
  private final File output;
  private final String ionMode;
  ExecuteMummiChog emc;
  

  public String getTaskDescription() {
    return "Calculating Mummichog... ";
  }

  public MummichogTask(PeakListRow[] peakListRows, double cutoff, String network,
      boolean force_primary_ion, String modeling,String ionMode, File output) {
    super();
    this.peakListRows = peakListRows;
    this.cutoff = cutoff;
    this.network = network;
    this.force_primary_ion = force_primary_ion;
    this.modeling = modeling;
    this.output = output;
    if(ionMode.equalsIgnoreCase("Positive-Default")) {
    	this.ionMode="pos_default";
    }else {
    	this.ionMode="negative";
    }
  }

  public double getFinishedPercentage() {
    if (emc != null) {
      return emc.getProgress() / 100 + finishedPercentage;
    } else {
      return finishedPercentage;
    }

  }

  public void run() {
    setStatus(TaskStatus.PROCESSING);
    String errorMsg = null;
    try {
      String input = prepareData();
      emc = new ExecuteMummiChog();

      MummichogParams mummiParameters = new MummichogParams(this.cutoff, this.network,
          this.force_primary_ion, this.modeling,this.ionMode, this.output.getAbsolutePath());


      System.out.println(input);
      @SuppressWarnings("unused")
      Map<String, List<Compound>> mummiOutput = emc.runMummiChog(input, mummiParameters);
      

      // Code to update the align peak list
      for (String mzr : mummiOutput.keySet()) {

        for (PeakListRow pr : peakListRows) {
          if ((String.valueOf(pr.getAverageMZ()) + ";" + String.valueOf(pr.getAverageRT()))
              .equalsIgnoreCase(mzr)) {
            for (Compound compund : mummiOutput.get(mzr)) {
              PeakIdentity pId =
                  new SimplePeakIdentity(compund.getName(), compund.getFormula(), "", "", "");
              pr.addPeakIdentity(pId, true);
            }
          }
        }
      }
      finishedPercentage = 0.90;
      setStatus(TaskStatus.FINISHED);
      logger.info("Calculating Mummichog is completed");
    } catch (Exception e) {
      errorMsg = "'Unknown Error' during Mummichog calculation: " + e.getMessage();
    } catch (Throwable t) {
      setStatus(TaskStatus.ERROR);
      setErrorMessage(t.getMessage());
      logger.log(Level.SEVERE, "Mummichog calculation error", t);
    }

    if (errorMsg != null) {
      setErrorMessage(errorMsg);
      setStatus(TaskStatus.ERROR);
    }
  }

  String prepareData() throws Exception {
    StringBuilder result = new StringBuilder();
    result.append(
        "m/z" + "\t" + "retention_time" + "\t" + "p-value" + "\t" + "t-score" + "\t" + "custom_id");
    result.append("\n");
    for (PeakListRow pr : this.peakListRows) {
      Map<String, String> properties = pr.getPeakInformation().getAllProperties();
      if (properties.containsKey(P_VALUE_KEY) && properties.containsKey(T_VALUE_KEY)) {
        result.append(pr.getAverageMZ()).append("\t").append(pr.getAverageRT()).append("\t")
            .append(pr.getPeakInformation().getAllProperties().get(P_VALUE_KEY)).append("\t")
            .append(pr.getPeakInformation().getAllProperties().get(T_VALUE_KEY)).append("\t")
            .append("randomText");
        result.append("\n");
      }
    }
    finishedPercentage = 0.10;
    return result.toString();
  }
}
