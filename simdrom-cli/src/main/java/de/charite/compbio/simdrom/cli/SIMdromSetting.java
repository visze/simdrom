package de.charite.compbio.simdrom.cli;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.MissingOptionException;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

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
	 * Probability so choose a mutation in the
	 * {@link SIMdromSetting#MUTATIONS_VCF}.
	 */
	public static double MUTATIONS_PROBABILITY = 1.0;
	/**
	 * If set only one random sample will be selected of the
	 * {@link SIMdromSetting#BACKGROUND_VCF}.
	 */
	public static boolean ONLY_ONE_SAMPLE = false;
	/**
	 * Identifier in the info-String of the allele frequency in the
	 * {@link SIMdromSetting#BACKGROUND_VCF} file.
	 */
	public static String ALLELE_FREQUENCY_IDENTIFIER;

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
		OptionBuilder
				.withDescription("Default 1.0. Choose variants with this probability. If larger than 1 the exact number of variants will be chosen (round down).");
		options.addOption(OptionBuilder.create());

		// mutations probability
		OptionBuilder.hasArg();
		OptionBuilder.withLongOpt("mutations-probability");
		OptionBuilder
				.withDescription("Default 1.0. Choose mutations with this probability. If larger than 1 the exact number of variants will be chosen (round down).");
		options.addOption(OptionBuilder.create());

		// only one sample
		OptionBuilder.withLongOpt("single-sample");
		OptionBuilder
				.withDescription("Default false. If present, a random sample will be chosen of the background VCF.");
		options.addOption(OptionBuilder.create());

		// only one sample
		OptionBuilder.hasArg();
		OptionBuilder.withLongOpt("allele-frequency-idenifier");
		OptionBuilder
				.withDescription("Optional. If set, the identifier in the info string of the background VCF will be used as single probabilities to call variants.");
		options.addOption(OptionBuilder.create("AF"));

		CommandLineParser parser = new GnuParser();
		try {
			CommandLine cmd = parser.parse(options, args);
			if (args.length == 0 || cmd.hasOption("h")) {
				throw new MissingOptionException("Please Insert an argument");
			}
			BACKGROUND_VCF = cmd.getOptionValue("background-population");
			if (cmd.hasOption("mutations"))
				MUTATIONS_VCF = cmd.getOptionValue("mutations");
			if (cmd.hasOption("background-probability"))
				BACKGROUND_PROBABILITY = Double.parseDouble(cmd.getOptionValue("background-probability"));
			if (cmd.hasOption("mutations-probability"))
				MUTATIONS_PROBABILITY = Double.parseDouble(cmd.getOptionValue("mutations-probability"));
			if (cmd.hasOption("single-sample"))
				ONLY_ONE_SAMPLE = true;
			if (cmd.hasOption("allele-frequency-idenifier"))
				ALLELE_FREQUENCY_IDENTIFIER = cmd.getOptionValue("allele-frequency-idenifier");

		} catch (MissingOptionException e) {
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("SIMdrom", options);
			System.exit(0);
		}
	}

}
