/*
 * Copyright 2006-2007 The MZmine Development Team
 * 
 * This file is part of MZmine.
 * 
 * MZmine is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 * 
 * MZmine is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * MZmine; if not, write to the Free Software Foundation, Inc., 51 Franklin St,
 * Fifth Floor, Boston, MA 02110-1301 USA
 */

package net.sf.mzmine.modules.alignment.join;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.util.logging.Logger;

import net.sf.mzmine.data.Parameter;
import net.sf.mzmine.data.ParameterSet;
import net.sf.mzmine.data.PeakList;
import net.sf.mzmine.data.Parameter.ParameterType;
import net.sf.mzmine.data.impl.SimpleParameter;
import net.sf.mzmine.data.impl.SimpleParameterSet;
import net.sf.mzmine.io.RawDataFile;
import net.sf.mzmine.main.MZmineCore;
import net.sf.mzmine.modules.BatchStep;
import net.sf.mzmine.project.MZmineProject;
import net.sf.mzmine.taskcontrol.Task;
import net.sf.mzmine.taskcontrol.TaskGroup;
import net.sf.mzmine.taskcontrol.TaskGroupListener;
import net.sf.mzmine.taskcontrol.TaskListener;
import net.sf.mzmine.userinterface.Desktop;
import net.sf.mzmine.userinterface.Desktop.MZmineMenu;
import net.sf.mzmine.userinterface.dialogs.ExitCode;
import net.sf.mzmine.userinterface.dialogs.ParameterSetupDialog;

/**
 * 
 */
public class JoinAligner implements BatchStep, TaskListener, ActionListener {

    public static final String RTToleranceTypeAbsolute = "Absolute";
    public static final String RTToleranceTypeRelative = "Relative";

    public static final Object[] RTToleranceTypePossibleValues = {
            RTToleranceTypeAbsolute, RTToleranceTypeRelative };

    public static final Parameter PeakListName = new SimpleParameter(
            ParameterType.STRING, "Peaklist name", "Peak list name");

    public static final Parameter MZvsRTBalance = new SimpleParameter(
            ParameterType.FLOAT, "M/Z vs RT balance",
            "Used in distance measuring as multiplier of M/Z difference", "",
            new Float(10.0), new Float(0.0), null);

    public static final Parameter MZTolerance = new SimpleParameter(
            ParameterType.FLOAT, "M/Z tolerance",
            "Maximum allowed M/Z difference", "Da", new Float(0.2), new Float(
                    0.0), null);

    public static final Parameter RTToleranceType = new SimpleParameter(
            ParameterType.STRING,
            "RT tolerance type",
            "Maximum RT difference can be defined either using absolute or relative value",
            RTToleranceTypeAbsolute, RTToleranceTypePossibleValues);

    public static final Parameter RTToleranceValueAbs = new SimpleParameter(
            ParameterType.FLOAT, "Absolute RT tolerance",
            "Maximum allowed absolute RT difference", "seconds",
            new Float(15.0), new Float(0.0), null);

    public static final Parameter RTToleranceValuePercent = new SimpleParameter(
            ParameterType.FLOAT, "Relative RT tolerance",
            "Maximum allowed relative RT difference", "%", new Float(0.15),
            new Float(0.0), null);

    private Logger logger = Logger.getLogger(this.getClass().getName());

    private ParameterSet parameters;

    private Desktop desktop;

    /**
     * @see net.sf.mzmine.main.MZmineModule#initModule(net.sf.mzmine.main.MZmineCore)
     */
    public void initModule() {

        this.desktop = MZmineCore.getDesktop();

        parameters = new SimpleParameterSet(new Parameter[] { PeakListName,
                MZvsRTBalance, MZTolerance, RTToleranceType,
                RTToleranceValueAbs, RTToleranceValuePercent });

        desktop.addMenuItem(MZmineMenu.ALIGNMENT, "Peak list aligner", this,
                null, KeyEvent.VK_A, false, true);

    }

    public String toString() {
        return new String("Join Aligner");
    }

    /**
     * @see net.sf.mzmine.main.MZmineModule#getParameterSet()
     */
    public ParameterSet getParameterSet() {
        return parameters;
    }

    public void setParameters(ParameterSet parameters) {
        this.parameters = parameters;
    }

    /**
     * @see net.sf.mzmine.modules.BatchStep#setupParameters(net.sf.mzmine.data.ParameterSet)
     */
    public ExitCode setupParameters(ParameterSet currentParameters) {
        ParameterSetupDialog dialog = new ParameterSetupDialog(
                desktop.getMainFrame(), "Please check parameter values for "
                        + toString(), (SimpleParameterSet) currentParameters);
        dialog.setVisible(true);
        return dialog.getExitCode();
    }

    /**
     * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
     */
    public void actionPerformed(ActionEvent e) {

        RawDataFile[] dataFiles = desktop.getSelectedDataFiles();
        MZmineProject currentProject = MZmineCore.getCurrentProject();

        if (dataFiles.length == 0) {
            desktop.displayErrorMessage("Please select at least one data file");
            return;
        }

        for (RawDataFile dataFile : dataFiles) {
            if (currentProject.getFilePeakList(dataFile) == null) {
                desktop.displayErrorMessage(dataFile
                        + " has no peak list. Please run peak picking first.");
                return;
            }
        }

        ExitCode exitCode = setupParameters(parameters);
        if (exitCode != ExitCode.OK)
            return;

        runModule(dataFiles, null, parameters.clone(), null);

    }

    public void taskStarted(Task task) {
        logger.info("Running join aligner");
    }

    public void taskFinished(Task task) {

        if (task.getStatus() == Task.TaskStatus.FINISHED) {

            logger.info("Finished join aligner");

            PeakList alignmentResult = (PeakList) task.getResult();

            MZmineCore.getCurrentProject().addAlignedPeakList(alignmentResult);

        } else if (task.getStatus() == Task.TaskStatus.ERROR) {
            /* Task encountered an error */
            String msg = "Error while aligning peak lists: "
                    + task.getErrorMessage();
            logger.severe(msg);
            desktop.displayErrorMessage(msg);

        }

    }

    /**
     * @see net.sf.mzmine.modules.BatchStep#runModule(net.sf.mzmine.io.RawDataFile[],
     *      net.sf.mzmine.data.PeakList[], net.sf.mzmine.data.ParameterSet,
     *      net.sf.mzmine.taskcontrol.TaskGroupListener)
     */
    public TaskGroup runModule(RawDataFile[] dataFiles,
            PeakList[] alignmentResults, ParameterSet parameters,
            TaskGroupListener methodListener) {

        MZmineProject currentProject = MZmineCore.getCurrentProject();

        // check peaklists
        for (int i = 0; i < dataFiles.length; i++) {
            if (currentProject.getFilePeakList(dataFiles[i]) == null) {
                String msg = "Cannot start alignment of " + dataFiles[i]
                        + ", please run peak picking first.";
                logger.severe(msg);
                desktop.displayErrorMessage(msg);
                return null;
            }
        }

        // prepare a new sequence with just one task
        Task tasks[] = new JoinAlignerTask[1];
        tasks[0] = new JoinAlignerTask(dataFiles,
                (SimpleParameterSet) parameters);
        TaskGroup newSequence = new TaskGroup(tasks, this, methodListener);

        // execute the sequence
        newSequence.run();

        return newSequence;

    }

}
