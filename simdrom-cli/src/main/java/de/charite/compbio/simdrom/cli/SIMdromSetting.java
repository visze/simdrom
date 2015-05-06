package de.charite.compbio.simdrom.cli;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.MissingOptionException;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import com.google.common.collect.ImmutableSet;

import de.charite.compbio.simdrom.filter.IFilter;
import de.charite.compbio.simdrom.filter.InfoFieldFilter;

public class SIMdromSetting {

	/**
	 * VCF of the background mutation. These mutations are used to sample a new
	 * mutation file. Required.
	 */
	public static String BACKGROUND_VCF;
	/**
	 * If set, mutations of these file are spiked in. Optional.
	 */
	public static String MUTATIONS_VCF;
	/**
	 * Probability so choose a variant in the
	 * {@link SIMdromSetting#BACKGROUND_VCF}.
	 */
	public static double BACKGROUND_PROBABILITY = 1.0;
	/**
	 * If set only the exact number of variants will be selected in the
	 * {@link SIMdromSetting#BACKGROUND_VCF}.
	 */
	public static int BACKGROUND_VARIANT_NUMBER;
	/**
	 * Probability so choose a mutation in the
	 * {@link SIMdromSetting#MUTATIONS_VCF}.
	 */
	public static double MUTATIONS_PROBABILITY = 1.0;
	/**
	 * If set only the exact number of variants will be selected in the
	 * {@link SIMdromSetting#MUTATIONS_VCF}.
	 */
	public static int MUTATIONS_VARIANT_NUMBER;
	/**
	 * If set only one random sample will be selected of the
	 * {@link SIMdromSetting#MUTATIONS_VCF}.
	 */
	public static boolean ONLY_ONE_SAMPLE = false;
	/**
	 * Identifier in the info-String of the allele frequency in the
	 * {@link SIMdromSetting#BACKGROUND_VCF} file.
	 */
	public static String BACKGROUND_ALLELE_FREQUENCY_IDENTIFIER;
	/**
	 * Identifier in the info-String of the allele frequency in the
	 * {@link SIMdromSetting#MUTATIONS_VCF} file.
	 */
	public static String MUTATIONS_ALLELE_FREQUENCY_IDENTIFIER;
	/**
	 * Identifier in the info-String of the ALT allele count in the
	 * {@link SIMdromSetting#MUTATIONS_VCF} file.
	 */
	public static String BACKGROUND_ALT_ALLELE_COUNT;
	/**
	 * Identifier in the info-String of the ALT allele count in the
	 * {@link SIMdromSetting#BACKGROUND_VCF} file.
	 */
	public static String MUTATIONS_ALT_ALLELE_COUNT;
	/**
	 * Identifier in the info-String of the all allele count in the
	 * {@link SIMdromSetting#MUTATIONS_VCF} file.
	 */
	public static String BACKGROUND_ALLELE_COUNT;
	/**
	 * Identifier in the info-String of the all allele count in the
	 * {@link SIMdromSetting#BACKGROUND_VCF} file.
	 */
	public static String MUTATIONS_ALLELE_COUNT;
	/**
	 * Spike in log file to get informations about the spike in.
	 */
	public static String SPLIKE_IN_LOGFILE;

	public static ImmutableSet<IFilter> MUTATIONS_FILTERS = ImmutableSet.<IFilter> builder().build();

	/**
	 * parse the option arguments of the command line and set the static fields.
	 * 
	 * @param args
	 * @throws ParseException
	 */
	public static void parse(String[] args) throws ParseException {
		Options options = new Options();

		// help
		OptionBuilder.withLongOpt("help");
		OptionBuilder.withDescription("Show this help message");
		options.addOption(OptionBuilder.create("h"));

		// background vcf
		OptionBuilder.hasArg();
		OptionBuilder.isRequired();
		OptionBuilder.withLongOpt("background-population");
		OptionBuilder
				.withDescription("VCF of the background population. variants will be uswed to sample a new mutation file.");
		options.addOption(OptionBuilder.create("b"));

		// mutations vcf
		OptionBuilder.hasArg();
		OptionBuilder.withLongOpt("mutations");
		OptionBuilder.withDescription("Optional. Mutation VCF to spike in.");
		options.addOption(OptionBuilder.create("m"));

		// background probability
		OptionBuilder.hasArg();
		OptionBuilder.withLongOpt("background-probability");
		OptionBuilder.withDescription("Default 1.0. Choose variants with this probability.");
		options.addOption(OptionBuilder.create());

		// background exact counts
		OptionBuilder.hasArg();
		OptionBuilder.withLongOpt("background-variants-amount");
		OptionBuilder
				.withDescription("Optional. Choose exact the given number of variants in the background population.");
		options.addOption(OptionBuilder.create());

		// mutations probability
		OptionBuilder.hasArg();
		OptionBuilder.withLongOpt("mutations-probability");
		OptionBuilder.withDescription("Default 1.0. Choose mutations with this probability.");
		options.addOption(OptionBuilder.create());

		// background exact counts
		OptionBuilder.hasArg();
		OptionBuilder.withLongOpt("mutations-variants-amount");
		OptionBuilder
				.withDescription("Optional. Choose exact the given number of variants in the mutation population.");
		options.addOption(OptionBuilder.create());

		// only one sample
		OptionBuilder.withLongOpt("single-sample");
		OptionBuilder
				.withDescription("Default false. If present, a random sample will be chosen of the background VCF.");
		options.addOption(OptionBuilder.create());

		// background allele frequency identifier
		OptionBuilder.hasArg();
		OptionBuilder.withLongOpt("background-allele-frequency-idenifier");
		OptionBuilder
				.withDescription("Optional. If set, the identifier in the info string of the background VCF will be used as single probabilities to call variants.");
		options.addOption(OptionBuilder.create("bAF"));

		// mutations allele frequency identifier
		OptionBuilder.hasArg();
		OptionBuilder.withLongOpt("mutations-allele-frequency-idenifier");
		OptionBuilder
				.withDescription("Optional. If set, the identifier in the info string of the mutations VCF will be used as single probabilities to call variants.");
		options.addOption(OptionBuilder.create("mAF"));

		// background allele ALT allele count
		OptionBuilder.hasArg();
		OptionBuilder.withLongOpt("background-alt-allele-count");
		OptionBuilder
				.withDescription("Optional. If set, the identifier in the info string of the background VCF will be used to compute single probabilities per variant. (bAC/bAN)");
		options.addOption(OptionBuilder.create("bAC"));

		// mutations allele ALT allele count
		OptionBuilder.hasArg();
		OptionBuilder.withLongOpt("mutations-alt-allele-count");
		OptionBuilder
				.withDescription("Optional. If set, the identifier in the info string of the mutations VCF will be used to compute single probabilities per variant. (mAC/mAN)");
		options.addOption(OptionBuilder.create("mAC"));
		// background allele allele count
		OptionBuilder.hasArg();
		OptionBuilder.withLongOpt("background-allele-count");
		OptionBuilder
				.withDescription("Optional. If set, the identifier in the info string of the background VCF will be used to compute single probabilities per variant. (bAC/bAN)");
		options.addOption(OptionBuilder.create("bAN"));

		// mutations allele allele count
		OptionBuilder.hasArg();
		OptionBuilder.withLongOpt("mutations-allele-count");
		OptionBuilder
				.withDescription("Optional. If set, the identifier in the info string of the mutations VCF will be used to compute single probabilities per variant. (mAC/mAN)");
		options.addOption(OptionBuilder.create("mAN"));

		// spike in log
		OptionBuilder.hasArg();
		OptionBuilder.withLongOpt("spike-in-log");
		OptionBuilder
				.withDescription("Optional. Path for a log file (TSV-Format) that descibes the spiked in mutations.");
		options.addOption(OptionBuilder.create());

		// mutations info filter
		OptionBuilder.hasArgs();
		OptionBuilder.withLongOpt("mutations-info-filter");
		OptionBuilder
				.withDescription("Optional. Uses the VCF info field to kepp only variants that passed the filter. Filter is written using the info field id followed by '=' and the value. Like CLNSIG=5");
		options.addOption(OptionBuilder.create());

		CommandLineParser parser = new GnuParser();
		try {
			CommandLine cmd = parser.parse(options, args);
			if (args.length == 0 || cmd.hasOption("h")) {
				throw new MissingOptionException("Please Insert an argument");
			}

			// check if input is correct
			checkNotAllowedOptions(cmd, "background-probability", "background-variants-amount",
					"background-allele-frequency-idenifier", "background-allele-count");
			checkNotAllowedOptions(cmd, "mutations-probability", "mutations-variants-amount",
					"mutations-allele-frequency-idenifier", "mutations-allele-count");
			checkMissingOption(cmd, "background-allele-count", "background-alt-allele-count");
			checkMissingOption(cmd, "mutations-allele-count", "mutations-alt-allele-count");

			BACKGROUND_VCF = cmd.getOptionValue("background-population");
			if (cmd.hasOption("mutations"))
				MUTATIONS_VCF = cmd.getOptionValue("mutations");

			// probabilities
			if (cmd.hasOption("background-probability")) {
				BACKGROUND_PROBABILITY = Double.parseDouble(cmd.getOptionValue("background-probability"));
			}
			if (cmd.hasOption("mutations-probability")) {
				MUTATIONS_PROBABILITY = Double.parseDouble(cmd.getOptionValue("mutations-probability"));
			}
			// variant counts
			if (cmd.hasOption("background-variants-amount")) {
				BACKGROUND_VARIANT_NUMBER = Integer.parseInt(cmd.getOptionValue("background-variants-amount"));
			}
			if (cmd.hasOption("mutations-variants-amount")) {
				MUTATIONS_VARIANT_NUMBER = Integer.parseInt(cmd.getOptionValue("mutations-variants-amount"));
			}
			// AF identifier
			if (cmd.hasOption("background-allele-frequency-idenifier")) {
				BACKGROUND_ALLELE_FREQUENCY_IDENTIFIER = cmd.getOptionValue("background-allele-frequency-idenifier");
			}
			if (cmd.hasOption("mutations-allele-frequency-idenifier")) {
				MUTATIONS_ALLELE_FREQUENCY_IDENTIFIER = cmd.getOptionValue("mutations-allele-frequency-idenifier");
			}
			// AC identifier
			if (cmd.hasOption("background-alt-allele-count")) {
				BACKGROUND_ALT_ALLELE_COUNT = cmd.getOptionValue("background-alt-allele-count");
			}
			if (cmd.hasOption("mutations-alt-allele-count")) {
				MUTATIONS_ALT_ALLELE_COUNT = cmd.getOptionValue("mutations-alt-allele-count");
			}
			// AN identifier
			if (cmd.hasOption("background-allele-count")) {
				BACKGROUND_ALLELE_COUNT = cmd.getOptionValue("background-allele-count");
			}
			if (cmd.hasOption("mutations-allele-count")) {
				MUTATIONS_ALLELE_COUNT = cmd.getOptionValue("mutations-allele-count");
			}
			// single sample
			if (cmd.hasOption("single-sample"))
				ONLY_ONE_SAMPLE = true;

			if (cmd.hasOption("spike-in-log"))
				SPLIKE_IN_LOGFILE = cmd.getOptionValue("spike-in-log");

			Set<IFilter> filters = new HashSet<IFilter>();
			if (cmd.hasOption("mutations-info-filter")) {
				for (String opt : cmd.getOptionValues("mutations-info-filter")) {
					String[] split = opt.split("=");
					if (isInt(split[1])) {
						filters.add(new InfoFieldFilter(split[0], Integer.parseInt(split[1])));
					} else if (isDouble(split[1])) {
						filters.add(new InfoFieldFilter(split[0], Double.parseDouble(split[1])));
					} else {
						filters.add(new InfoFieldFilter(split[0], split[1]));
					}
				}
			}

			MUTATIONS_FILTERS = ImmutableSet.<IFilter> builder().addAll(filters).build();
		} catch (MissingOptionException e) {
			HelpFormatter formatter = new HelpFormatter();
			formatter.setWidth(120);
			formatter.printHelp("SIMdrom", options);
			System.exit(0);
		} catch (NotAllowedCombinationOfOptionsException | MissingOptionsException e) {
			e.printStackTrace();
			System.exit(0);
		}
	}

	private static void checkNotAllowedOptions(CommandLine cmd, String... values)
			throws NotAllowedCombinationOfOptionsException {
		List<String> falseOptions = new ArrayList<String>();
		for (String opt : values) {
			if (cmd.hasOption(opt))
				falseOptions.add(opt);
		}
		if (falseOptions.size() > 1)
			throw new NotAllowedCombinationOfOptionsException(falseOptions);
	}

	private static void checkMissingOption(CommandLine cmd, String... values) throws MissingOptionsException {
		List<String> missingOptions = new ArrayList<String>();
		List<String> usedOptions = new ArrayList<String>();
		for (String opt : values) {
			if (cmd.hasOption(opt))
				usedOptions.add(opt);
			else
				missingOptions.add(opt);
		}
		if (usedOptions.size() > 0 && missingOptions.size() > 0)
			throw new MissingOptionsException(usedOptions, missingOptions);
	}

	private static boolean isDouble(String string) {
		return Pattern.matches("^\\d+.\\d+$", string);
	}

	private static boolean isInt(String string) {
		return Pattern.matches("^\\d+$", string);
	}

}
