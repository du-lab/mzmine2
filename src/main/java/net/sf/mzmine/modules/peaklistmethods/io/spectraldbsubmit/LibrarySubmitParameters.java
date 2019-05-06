/*
 * Copyright 2006-2018 The MZmine 2 Development Team
 * 
 * This file is part of MZmine 2.
 * 
 * MZmine 2 is free software; you can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 * 
 * MZmine 2 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
 * even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with MZmine 2; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
 * USA
 */

/*
 * This module was prepared by Abi Sarvepalli, Christopher Jensen, and Zheng Zhang at the Dorrestein
 * Lab (University of California, San Diego).
 * 
 * It is freely available under the GNU GPL licence of MZmine2.
 * 
 * For any questions or concerns, please refer to:
 * https://groups.google.com/forum/#!forum/molecular_networking_bug_reports
 * 
 * Credit to the Du-Lab development team for the initial commitment to the MGF export module.
 */

package net.sf.mzmine.modules.peaklistmethods.io.spectraldbsubmit;

import java.text.DecimalFormat;
import net.sf.mzmine.main.MZmineCore;
import net.sf.mzmine.modules.peaklistmethods.io.spectraldbsubmit.GnpsValues.CompoundSource;
import net.sf.mzmine.modules.peaklistmethods.io.spectraldbsubmit.GnpsValues.Instrument;
import net.sf.mzmine.modules.peaklistmethods.io.spectraldbsubmit.GnpsValues.IonSource;
import net.sf.mzmine.modules.peaklistmethods.io.spectraldbsubmit.GnpsValues.Polarity;
import net.sf.mzmine.parameters.Parameter;
import net.sf.mzmine.parameters.impl.SimpleParameterSet;
import net.sf.mzmine.parameters.parametertypes.BooleanParameter;
import net.sf.mzmine.parameters.parametertypes.ComboParameter;
import net.sf.mzmine.parameters.parametertypes.DoubleParameter;
import net.sf.mzmine.parameters.parametertypes.OptionalParameter;
import net.sf.mzmine.parameters.parametertypes.PasswordParameter;
import net.sf.mzmine.parameters.parametertypes.StringParameter;
import net.sf.mzmine.parameters.parametertypes.filenames.FileNameParameter;

/**
 * 
 * @author Robin Schmid (robinschmid@uni-muenster.de)
 *
 */
public class LibrarySubmitParameters extends SimpleParameterSet {

  /**
   * TODO difference of COMPOUND MASS and EXACT mass
   * 
   * @author r_schm33
   *
   */


  public static final ComboParameter<CompoundSource> ACQUISITION =
      new ComboParameter<>("ACQUISITION", "", CompoundSource.values(), CompoundSource.Crude);
  public static final ComboParameter<Polarity> IONMODE =
      new ComboParameter<>("IONMODE", "", Polarity.values(), Polarity.Positive);
  public static final ComboParameter<Instrument> INSTRUMENT =
      new ComboParameter<>("INSTRUMENT", "", Instrument.values(), Instrument.Orbitrap);
  public static final ComboParameter<IonSource> ION_SOURCE =
      new ComboParameter<>("IONSOURCE", "", IonSource.values(), IonSource.LC_ESI);

  // save to local file
  public static final OptionalParameter<FileNameParameter> LOCALFILE =
      new OptionalParameter<>(new FileNameParameter("Local file", "Local library file"), false);
  public static final BooleanParameter EXPORT_GNPS_JSON = new BooleanParameter(
      "Export GNPS json file", "The GNPS library submission json format", true);
  public static final BooleanParameter EXPORT_MSP =
      new BooleanParameter("Export NIST msp file", "The NIST msp library format", true);
  // user and password
  public static final BooleanParameter SUBMIT_GNPS =
      new BooleanParameter("Submit to GNPS", "Submit new entry to GNPS library", true);
  public static final StringParameter USERNAME = new StringParameter("username", "", "", true);
  public static final PasswordParameter PASSWORD = new PasswordParameter("password", "", "", true);
  // all general parameters
  public static final StringParameter DESCRIPTION =
      new StringParameter("description", "", "", false);
  public static final StringParameter COMPOUND_NAME =
      new StringParameter("COMPOUND_NAME", "", "", true);
  public static final StringParameter PI =
      new StringParameter("PI", "Principle investigator", "", true);
  public static final StringParameter DATACOLLECTOR =
      new StringParameter("DATACOLLECTOR", "", "", true);

  public static final StringParameter FRAGMENTATION_METHOD =
      new StringParameter("FRAGMENTATION_METHOD", "", "", false);
  public static final StringParameter INSTRUMENT_NAME =
      new StringParameter("INSTRUMENT_NAME", "", "", false);
  public static final StringParameter DATA_COLLECTOR =
      new StringParameter("DATACOLLECTOR", "", "", false);
  public static final StringParameter PUBMED = new StringParameter("PUBMED", "", "", false);
  public static final StringParameter INCHI_AUX = new StringParameter("INCHIAUX", "", "", false);
  public static final StringParameter INCHI =
      new StringParameter("INCHI", "Structure as INCHI", "", false);
  public static final StringParameter CAS = new StringParameter("CASNUMBER", "", "", false);
  public static final StringParameter SMILES =
      new StringParameter("SMILES", "Structure as SMILES code", "", false);
  public static final StringParameter FORMULA = new StringParameter("FORMULA", "", "", false);
  public static final DoubleParameter EXACT_MASS = new DoubleParameter("EXACTMASS",
      "Monoisotopic neutral mass of compound", new DecimalFormat("0.000###"), 0d);


  public static final OptionalParameter<DoubleParameter> EXPORT_RT = new OptionalParameter<>(
      new DoubleParameter("RT", "Retention time", MZmineCore.getConfiguration().getRTFormat(), 0d),
      false);

  // is not used: this would override MZ (the precursor MZ)
  // public static final DoubleParameter MOLECULE_MASS = new DoubleParameter("MOLECULEMASS",
  // "Exact precursor m/z", MZmineCore.getConfiguration().getMZFormat(), 0d);

  public LibrarySubmitParameters() {
    super(new Parameter[] {
        // save to local file
        LOCALFILE, EXPORT_GNPS_JSON, EXPORT_MSP,
        // submit to online library
        SUBMIT_GNPS,
        // username password
        USERNAME, PASSWORD,
        // Always set
        PI, DATA_COLLECTOR, DESCRIPTION, COMPOUND_NAME, EXACT_MASS, EXPORT_RT, INSTRUMENT_NAME,
        INSTRUMENT, ION_SOURCE, ACQUISITION, IONMODE, FRAGMENTATION_METHOD,
        // newly introduced
        FORMULA,
        // optional
        SMILES, INCHI, INCHI_AUX, CAS, PUBMED});
  }
}
