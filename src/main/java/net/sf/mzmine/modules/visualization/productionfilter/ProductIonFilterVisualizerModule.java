package net.sf.mzmine.modules.visualization.productionfilter;

import java.util.Collection;

import javax.annotation.Nonnull;

import net.sf.mzmine.datamodel.MZmineProject;
import net.sf.mzmine.datamodel.RawDataFile;
import net.sf.mzmine.modules.MZmineModuleCategory;
import net.sf.mzmine.modules.MZmineRunnableModule;
import net.sf.mzmine.parameters.ParameterSet;
import net.sf.mzmine.taskcontrol.Task;
import net.sf.mzmine.util.ExitCode;

/**
 * Production ion filter visualizer generated by Shawn Hoogstra : shoogstr@uwo.ca
 * Filter ms/ms scan based upon desired m/z and neutral loss m/z values
 */
public class ProductIonFilterVisualizerModule implements MZmineRunnableModule {

	private static final String MODULE_NAME = "Diagnostic fragmentation filtering";
	private static final String MODULE_DESCRIPTION = "This visualizer will filter MS/MS scans based on desired m/z and neutral loss and output a"
			+ "plot showing the product ions vs. the precursor ion. In addition it will output a file containing m/z and RT of filtered MS/MS scans. ";

	@Override
	public @Nonnull String getName() {
		return MODULE_NAME;
	}

	@Override
	public @Nonnull String getDescription() {
		return MODULE_DESCRIPTION;
	}

	@Override
	@Nonnull
	public ExitCode runModule(@Nonnull MZmineProject project, @Nonnull ParameterSet parameters,
			@Nonnull Collection<Task> tasks) {

		RawDataFile dataFiles[] = parameters.getParameter(ProductIonFilterParameters.dataFiles).getValue()
				.getMatchingRawDataFiles();
		
		
		ProductIonFilterVisualizerWindow newWindow = new ProductIonFilterVisualizerWindow(dataFiles[0], parameters);
		newWindow.setVisible(true);

		return ExitCode.OK;
	}

	@Override
	public @Nonnull MZmineModuleCategory getModuleCategory() {
		return MZmineModuleCategory.VISUALIZATIONRAWDATA;
	}

	@Override
	public @Nonnull Class<? extends ParameterSet> getParameterSetClass() {
		return ProductIonFilterParameters.class;
	}

}