/*
 * Copyright (C) 2018 Du-Lab Team <dulab.binf@gmail.com>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
package net.sf.mzmine.modules.peaklistmethods.peakpicking.adapsimplespectraldeconvolution;

import net.sf.mzmine.datamodel.MZmineProject;
import net.sf.mzmine.datamodel.PeakList;
import net.sf.mzmine.datamodel.RawDataFile;
import net.sf.mzmine.modules.MZmineModuleCategory;
import net.sf.mzmine.modules.MZmineProcessingModule;
import net.sf.mzmine.modules.peaklistmethods.peakpicking.adap3decompositionV2.ChromatogramPeakPair;
import net.sf.mzmine.parameters.ParameterSet;
import net.sf.mzmine.taskcontrol.Task;
import net.sf.mzmine.util.ExitCode;

import javax.annotation.Nonnull;
import java.util.Collection;
import java.util.Map;

/**
 *
 * @author aleksandrsmirnov
 */
public class SimpleSpectralDeconvolutionModule implements MZmineProcessingModule {
    
    private static final String MODULE_NAME = "Simple Spectral Deconvolution";
    private static final String MODULE_DESCRIPTION = "This method "
            + "combines peaks with similar EICs";

    @Override
    public @Nonnull String getName() {
        return MODULE_NAME;
    }

    @Override
    public @Nonnull String getDescription() {
        return MODULE_DESCRIPTION;
    }
    
    @Override
    public @Nonnull MZmineModuleCategory getModuleCategory() {
        return MZmineModuleCategory.SPECTRALDECONVOLUTION;
    }

    @Override
    public @Nonnull Class<? extends ParameterSet> getParameterSetClass() {
        return SimpleSpectralDeconvolutionParameters.class;
    }
    
    @Override @Nonnull
    public ExitCode runModule(@Nonnull MZmineProject project,
            @Nonnull ParameterSet parameters, @Nonnull Collection<Task> tasks)
    {
        PeakList[] peakLists = parameters
                .getParameter(SimpleSpectralDeconvolutionParameters.PEAK_LISTS).getValue()
                .getMatchingPeakLists();

        for (PeakList peakList : peakLists) {
            Task newTask = new SimpleSpectralDeconvolutionTask(project, peakList, parameters);
            tasks.add(newTask);
        }
        
        return ExitCode.OK;
    }
}