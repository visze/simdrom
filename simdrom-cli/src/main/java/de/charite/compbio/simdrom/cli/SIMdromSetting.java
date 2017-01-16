package de.charite.compbio.simdrom.cli;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.MissingOptionException;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import com.google.common.collect.ImmutableSet;

import de.charite.compbio.simdrom.cli.exception.MissingOptionsException;
import de.charite.compbio.simdrom.cli.exception.NotAllowedCombinationOfOptionsException;
import de.charite.compbio.simdrom.cli.exception.WrongIntervalFormatException;
import de.charite.compbio.simdrom.filter.EqualInfoFieldFilter;
import de.charite.compbio.simdrom.filter.IFilter;
import de.charite.compbio.simdrom.interval.SAMFileHeaderBuilder;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;

/**
 * Command line options class for the SIMdrom.
 * 
 * @author <a href="mailto:max.schubach@charite.de">Max Schubach</a>
 *
 */
public class SIMdromSetting {

	/**
	 * VCF of the background mutation. These mutations are used to sample a new mutation file. Required.
	 */
	public static File BACKGROUND_VCF;
	/**
	 * If set, mutations of these file are spiked in. Optional.
	 */
	public static File MUTATIONS_VCF;
	/**
	 * Probability so choose a variant in the {@link SIMdromSetting#BACKGROUND_VCF}.
	 */
	public static double BACKGROUND_PROBABILITY = 1.0;
	/**
	 * If set only the exact number of variants will be selected in the {@link SIMdromSetting#BACKGROUND_VCF}.
	 */
	public static int BACKGROUND_VARIANT_NUMBER;
	/**
	 * Probability so choose a mutation in the {@link SIMdromSetting#MUTATIONS_VCF}.
	 */
	public static double MUTATIONS_PROBABILITY = 1.0;
	/**
	 * If set only the exact number of variants will be selected in the {@link SIMdromSetting#MUTATIONS_VCF}.
	 */
	public static int MUTATIONS_VARIANT_NUMBER;
	/**
	 * If set only one (random) sample will be selected of the {@link SIMdromSetting#BACKGROUND_VCF}.
	 */
	public static boolean ONLY_ONE_BACKGROUND_SAMPLE = false;
	/**
	 * If not null the given sample will be selected {@link SIMdromSetting#BACKGROUND_VCF}.
	 */
	public static String ONLY_ONE_BACKGROUND_SAMPLE_NAME;
	/**
	 * If set only one (random) sample will be selected of the {@link SIMdromSetting#MUTATIONS_VCF}.
	 */
	public static boolean ONLY_ONE_MUTATIONS_SAMPLE = false;
	/**
	 * If not null the given sample will be selected {@link SIMdromSetting#MUTATIONS_VCF}.
	 */
	public static String ONLY_ONE_MUTATIONS_SAMPLE_NAME;
	/**
	 * Identifier in the info-String of the allele frequency in the {@link SIMdromSetting#BACKGROUND_VCF} file.
	 */
	public static String BACKGROUND_ALLELE_FREQUENCY_IDENTIFIER;
	/**
	 * Identifier in the info-String of the allele frequency in the {@link SIMdromSetting#MUTATIONS_VCF} file.
	 */
	public static String MUTATIONS_ALLELE_FREQUENCY_IDENTIFIER;
	/**
	 * Identifier in the info-String of the ALT allele count in the {@link SIMdromSetting#MUTATIONS_VCF} file.
	 */
	public static String BACKGROUND_ALT_ALLELE_COUNT;
	/**
	 * Identifier in the info-String of the ALT allele count in the {@link SIMdromSetting#BACKGROUND_VCF} file.
	 */
	public static String MUTATIONS_ALT_ALLELE_COUNT;
	/**
	 * Identifier in the info-String of the all allele count in the {@link SIMdromSetting#MUTATIONS_VCF} file.
	 */
	public static String BACKGROUND_ALLELE_COUNT;
	/**
	 * Identifier in the info-String of the all allele count in the {@link SIMdromSetting#BACKGROUND_VCF} file.
	 */
	public static String MUTATIONS_ALLELE_COUNT;
	/**
	 * If set, generate deNovo mutations.
	 */
	public static boolean USE_DE_NOVO = false;
	/**
	 * DeNovo rate or new mutations.
	 */
	public static double DE_NOVO_RATE = 1.2 * Math.pow(10, -8);
	/**
	 * Reference file
	 */
	public static String REFERENCE;
	/**
	 * Spike in log file to get informations about the spike in.
	 */
	public static String SPLIKE_IN_LOGFILE;
	/**
	 * Intervals. only write out at these points.
	 */
	public static Optional<IntervalList> INTERVALS;
	/**
	 * Output file. null if standard out.
	 */
	public static String OUTPUT;
	/**
	 * Mutation filter
	 */
	public static ImmutableSet<IFilter> MUTATIONS_FILTERS;

	/**
	 * parse the option arguments of the command line and set the static fields.
	 * 
	 * @param args
	 *            Arguments of options
	 * @throws ParseException
	 *             if arguments coud not parse
	 */
	public static void parse(String[] args) throws ParseException {
		Options options = new Options();

		// help
		options.addOption(Option.builder("h").longOpt("help").desc("Show this help message").build());

		// background vcf
		options.addOption(Option.builder("b").longOpt("background-population").hasArg().required().type(File.class)
				.desc("VCF of the background population. variants will be used to sample a new mutation file.")
				.build());

		// mutations vcf
		options.addOption(Option.builder("m").longOpt("mutations").hasArg().type(File.class)
				.desc("Optional. Mutation VCF to spike in.").build());

		// background probability
		options.addOption(Option.builder().longOpt("background-probability").hasArg()
				.desc("Default 1.0. Choose variants with this probability.").build());

		// background exact counts
		options.addOption(Option.builder().longOpt("background-variants-amount").hasArg()
				.desc("Optional. Choose exact the given number of variants in the background population.").build());

		// mutations probability
		options.addOption(Option.builder().longOpt("mutations-probability").hasArg()
				.desc("Default 1.0. Choose mutations with this probability.").build());

		// background exact counts
		options.addOption(Option.builder().longOpt("mutations-variants-amount").hasArg()
				.desc("Optional. Choose exact the given number of variants in the mutation population.").build());

		// only one sample of background
		options.addOption(Option.builder().longOpt("single-sample-background").hasArg().optionalArg(true)
				.desc("Default false. If present, a random sample will be chosen of the background VCF (-b).").build());

		// only one sample of mutations
		options.addOption(Option.builder().longOpt("single-sample-mutations").hasArg().optionalArg(true)
				.desc("Default false. If present, a random sample will be chosen of the mutations VCF (-m).").build());

		// background allele frequency identifier
		options.addOption(Option.builder("bAF").longOpt("background-allele-frequency-identifier").hasArg()
				.desc("Optional. If set, the identifier in the info string of the background VCF will be used as single probabilities to call variants.")
				.build());

		// mutations allele frequency identifier
		options.addOption(Option.builder("mAF").hasArg().longOpt("mutations-allele-frequency-identifier")
				.desc("Optional. If set, the identifier in the info string of the mutations VCF will be used as single probabilities to call variants.")
				.build());

		// background allele ALT allele count
		options.addOption(Option.builder("bAC").hasArg().longOpt("background-alt-allele-count")
				.desc("Optional. If set, the identifier in the info string of the background VCF will be used to compute single probabilities per variant. (bAC/bAN)")
				.build());

		// mutations allele ALT allele count
		options.addOption(Option.builder("mAC").hasArg().longOpt("mutations-alt-allele-count")
				.desc("Optional. If set, the identifier in the info string of the mutations VCF will be used to compute single probabilities per variant. (mAC/mAN)")
				.build());

		// background allele allele count
		options.addOption(Option.builder("bAN").hasArg().longOpt("background-allele-count")
				.desc("Optional. If set, the identifier in the info string of the background VCF will be used to compute single probabilities per variant. (bAC/bAN)")
				.build());

		// mutations allele allele count
		options.addOption(Option.builder("mAN").hasArg().longOpt("mutations-allele-count")
				.desc("Optional. If set, the identifier in the info string of the mutations VCF will be used to compute single probabilities per variant. (mAC/mAN)")
				.build());

		// deNovo rate
		options.addOption(Option.builder().optionalArg(true).longOpt("de-novo")
				.desc("Optional. If set, de-novo mutations are spiked in. Standard rate is 1.2*10^-8. But you can provide your own rate with this option. An indexed reference have to be set (see option --reference).")
				.build());

		// Reference file
		options.addOption(Option.builder().hasArg().longOpt("reference")
				.desc("Needed for option --de-novo. Please enter the paths to an indexed multi-FASTA file of your reference genome.")
				.build());

		// spike in log
		options.addOption(Option.builder().hasArg().longOpt("spike-in-log")
				.desc("Optional. Path for a log file (TSV-Format) that descibes the spiked in mutations.").build());

		// spike in log
		options.addOption(Option.builder("i").hasArgs().longOpt("interval")
				.desc("Optional. Use the parameter with intervals (chr1:12113-12123) or insert an interval list (file constist of one interval in each line)")
				.build());

		// mutations info filter
		options.addOption(Option.builder().hasArgs().longOpt("mutations-info-filter")
				.desc("Optional. Uses the VCF info field to kepp only variants that passed the filter. Filter is written using the info field id followed by '=' and the value. Like CLNSIG=5")
				.build());

		// output
		options.addOption(Option.builder("o").hasArg().longOpt("output")
				.desc("Optional. Writes the variants into this (bgzip) VCF file instead of printing it to the standard output.")
				.build());

		CommandLineParser parser = new DefaultParser();
		try {
			CommandLine cmd = parser.parse(options, args);
			if (args.length == 0 || cmd.hasOption("h")) {
				throw new MissingOptionException("Please Insert an argument");
			}

			// check if input is correct
			checkNotAllowedOptions(cmd, "background-probability", "background-variants-amount",
					"background-allele-frequency-identifier", "background-allele-count");
			checkNotAllowedOptions(cmd, "mutations-probability", "mutations-variants-amount",
					"mutations-allele-frequency-identifier", "mutations-allele-count");
			checkMissingOption(cmd, "background-allele-count", "background-alt-allele-count");
			checkMissingOption(cmd, "mutations-allele-count", "mutations-alt-allele-count");
			checkMissingOption(cmd, "de-novo", "reference");

			BACKGROUND_VCF = (File) cmd.getParsedOptionValue("background-population");
			if (cmd.hasOption("mutations"))
				MUTATIONS_VCF = (File) cmd.getParsedOptionValue("mutations");

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
			if (cmd.hasOption("background-allele-frequency-identifier")) {
				BACKGROUND_ALLELE_FREQUENCY_IDENTIFIER = cmd.getOptionValue("background-allele-frequency-identifier");
			}
			if (cmd.hasOption("mutations-allele-frequency-identifier")) {
				MUTATIONS_ALLELE_FREQUENCY_IDENTIFIER = cmd.getOptionValue("mutations-allele-frequency-identifier");
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
			// single sample background
			if (cmd.hasOption("single-sample-background")) {
				ONLY_ONE_BACKGROUND_SAMPLE = true;
				ONLY_ONE_BACKGROUND_SAMPLE_NAME = cmd.getOptionValue("single-sample-background");
			}
			// single sample mutations
			if (cmd.hasOption("single-sample-mutations")) {
				ONLY_ONE_MUTATIONS_SAMPLE = true;
				ONLY_ONE_MUTATIONS_SAMPLE_NAME = cmd.getOptionValue("single-sample-mutations");
			}
			// de novo
			if (cmd.hasOption("de-novo")) {
				USE_DE_NOVO = true;
				if (cmd.getOptionValue("de-novo") != null) {
					DE_NOVO_RATE = Double.parseDouble(cmd.getOptionValue("de-novo"));
				}
				REFERENCE = cmd.getOptionValue("reference");
			}
			// spike in log
			if (cmd.hasOption("spike-in-log"))
				SPLIKE_IN_LOGFILE = cmd.getOptionValue("spike-in-log");
			// intervals
			if (cmd.hasOption("interval")) {
				List<Interval> lst = new ArrayList<Interval>();
				for (String intervalString : cmd.getOptionValues("interval")) {
					lst.addAll(getIntervalOfOption(intervalString));
				}
				IntervalList list = new IntervalList(SAMFileHeaderBuilder.build());
				list.addall(lst);
				INTERVALS = Optional.of(list);
			} else
				INTERVALS = Optional.empty();
			// filters
			Set<IFilter> filters = new HashSet<IFilter>();
			if (cmd.hasOption("mutations-info-filter")) {
				for (String opt : cmd.getOptionValues("mutations-info-filter")) {
					String[] split = opt.split("=");
					if (isInt(split[1])) {
						filters.add(new EqualInfoFieldFilter(split[0], Integer.parseInt(split[1])));
					} else if (isDouble(split[1])) {
						filters.add(new EqualInfoFieldFilter(split[0], Double.parseDouble(split[1])));
					} else {
						filters.add(new EqualInfoFieldFilter(split[0], split[1]));
					}
				}
			}
			MUTATIONS_FILTERS = ImmutableSet.<IFilter> builder().addAll(filters).build();

			// output
			if (cmd.hasOption("output")) {
				OUTPUT = cmd.getOptionValue("output");
			}
		} catch (MissingOptionException e) {
			HelpFormatter formatter = new HelpFormatter();
			formatter.setWidth(120);
			formatter.printHelp("SIMdrom", options);
			System.exit(0);
		} catch (NotAllowedCombinationOfOptionsException | MissingOptionsException | IOException
				| WrongIntervalFormatException e) {
			e.printStackTrace();
			System.exit(0);
		}
	}

	private static List<Interval> getIntervalOfOption(String intervalString)
			throws IOException, WrongIntervalFormatException {
		try {
			List<Interval> output = new ArrayList<Interval>();
			output.add(getInterval(intervalString));
			return output;
		} catch (WrongIntervalFormatException e) {
			return getIntervalsOfFile(intervalString);
		}

	}

	private static Interval getInterval(String intervalString) throws WrongIntervalFormatException {
		Pattern intervalPattern = Pattern.compile("^((chr)?(\\d+|[XYM])):(\\d+)-(\\d+)$");
		Matcher m = intervalPattern.matcher(intervalString);
		if (m.matches()) {
			return new Interval(m.group(1), Integer.parseInt(m.group(4)), Integer.parseInt(m.group(5)));
		}
		throw new WrongIntervalFormatException(intervalString);
	}

	private static List<Interval> getIntervalsOfFile(String filepath) throws IOException, WrongIntervalFormatException {
		List<Interval> output = new ArrayList<Interval>();
		File file = new File(filepath);
		BufferedReader br = new BufferedReader(new FileReader(file));
		String line;
		while ((line = br.readLine()) != null) {
			if (line.trim().isEmpty())
				continue;
			output.add(getInterval(line.trim()));
		}
		br.close();
		return output;
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
