/**
 * 
 */
package de.charite.compbio.simdrom.sampler.vcf;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;
import java.util.Random;
import java.util.Set;
import java.util.stream.Stream;

import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Lists;
import com.google.common.collect.MoreCollectors;

import de.charite.compbio.simdrom.exception.InfoIDNotFoundException;
import de.charite.compbio.simdrom.filter.IFilter;
import de.charite.compbio.simdrom.sampler.DeNovoSampler;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;

/**
 * 
 * The {@link VCFSampler} class uses a VCF file as basis and acts as an iterator over such a file. There are multiple
 * possible options what the {@link VCFSampler} does with the read {@link VariantContext}s.
 * <p>
 * <dl>
 * <dt>Filter</dt>
 * <dd>Variants can be filtered using multiple {@link IFilter}. Completely removed variants where not finally passed
 * into the iterator.</dd>
 * <dt>Sample selection</dt>
 * <dd>If a multi-VCF file is used there is the possibility to select one sample out of it. It can be done randomly by
 * running the VCFRandomSampleSelecter or you can predefine a sample name. The sample name must be present in the VCF
 * header.</dd>
 * <dt>Variant generation using Hardy-Weinberg principle</dt>
 * <dd>If allele frequencies are present in the pares VCF they can be used to generate a completely new artificial
 * {@link Genotype} by using the Hardy-Weinberg principle. Direct allele frequencies can be used of the allele count
 * divided by the number of total alleles.</dd> *
 * <dt>Interval list</dt>
 * <dd>If an {@link IntervalList} is present only variants of that lists will be pares by the {@link VCFFileReader} and
 * these variants are only used by the {@link VCFSampler}.</dd>
 * <dt>De-novo generation</dt>
 * <dd>Inventing new variations. This function is not implemented yet.</dd>
 * </dl>
 * 
 * @author <a href="mailto:max.schubach@charite.de">Max Schubach</a>
 *
 */
public class VCFSampler implements CloseableIterator<VariantContext> {

	/**
	 * The default sample name for genotypes in the VCF.
	 */
	public static final String DEFAULT_SAMPLE_NAME = "Sampled";
	/**
	 * The probability to select a variant. If set to 1 all variants will be selected.
	 */
	private double probability;
	/**
	 * 
	 */
	private List<Integer> selectAlleles;
	/**
	 * The maximum number of variants that should be selected by the sampler. If set to smaller than 1 this option is of
	 * and there will be no maximum limit used.
	 */
	private int variantsAmount;
	private int position = -1;
	/**
	 * The identifier of the allele frequency in the Info column.
	 */
	private String afIdentifier;
	/**
	 * The identifier of the alt allele counts in the Info column.
	 */
	private String acIdentifier;
	/**
	 * The identifier of the total number of alleles.
	 */
	private String anIdentifier;
	/**
	 * If set the genotypes of this samples are selected.
	 */
	private Optional<String> sample;
	/**
	 * The VCF file reader.
	 */
	private VCFFileReader reader;

	/**
	 * The variant iterator defined by the reader and maybe by intervals.
	 */
	private CloseableIterator<VariantContext> iterator;
	/**
	 * The next variantContext. If this is null there will be no further variant and {@link #hasNext()} is false.
	 */
	private VariantContext next;
	/**
	 * Random number generator. Used if probabilities to choose a variant are smaller than 1.0
	 */
	private Random random;
	/**
	 * De novo generator TODO not used
	 */
	@SuppressWarnings("unused")
	private DeNovoSampler deNovoGenerator;
	/**
	 * Set of filters that will be used for each {@link VariantContext} read by the {@link #reader}
	 */
	private ImmutableSet<IFilter> filters;
	/**
	 * If list is set (not null) for each interval the {@link #reader} is queried using the intervals and only such
	 * variants are parsed
	 */
	private IntervalList intervals;
	/**
	 * Points to the actual interval in the interval list. If larger than {@link IntervalList#size()} then variants from
	 * all intervals are parsed.
	 */
	private int intervalPosition = 0;

	/**
	 * Default constructor initializing all fields. Can only be used by the {@link Builder}.
	 * 
	 * @param reader
	 *            {@link #probability}
	 * @param probability
	 *            {@link #probability}
	 * @param sample
	 *            {@link #sample}
	 * @param variantsAmount
	 *            {@link #variantsAmount}
	 * @param afIdentifier
	 *            {@link #afIdentifier}
	 * @param acIdentifier
	 *            {@link #acIdentifier}
	 * @param anIdentifier
	 *            {@link #anIdentifier}
	 * @param filters
	 *            {@link #filters}
	 * @param intervals
	 *            {@link #intervals}
	 * @param deNovoSampler
	 *            {@link #deNovoGenerator}
	 * @param selectAlleles
	 *            {@link #selectAlleles}
	 * @param random
	 *            {@link #random}
	 */
	private VCFSampler(VCFFileReader reader, double probability, Optional<String> sample, int variantsAmount,
			String afIdentifier, String acIdentifier, String anIdentifier, ImmutableSet<IFilter> filters,
			IntervalList intervals, DeNovoSampler deNovoSampler, List<Integer> selectAlleles, Random random) {
		this.reader = reader;
		this.probability = probability;
		this.sample = sample;
		this.variantsAmount = variantsAmount;
		this.afIdentifier = afIdentifier;
		this.acIdentifier = acIdentifier;
		this.anIdentifier = anIdentifier;
		this.filters = filters;
		this.intervals = intervals;
		this.selectAlleles = selectAlleles;
		this.deNovoGenerator = deNovoSampler;
		this.random = random;

		next();
	}

	/**
	 * 
	 * Builder class of the {@link VCFSampler}. Checks if everything is set correctly before building. The reader have
	 * to be set by using {@link #vcfReader(VCFFileReader)}.
	 * 
	 * @author <a href="mailto:max.schubach@charite.de">Max Schubach</a>
	 *
	 */
	public static final class Builder {

		/**
		 * See {@link VCFSampler#reader}
		 */
		private VCFFileReader reader;
		/**
		 * See {@link VCFSampler#probability}
		 */
		private double probability = 1.0;
		private List<Integer> selectAlleles;
		/**
		 * See {@link VCFSampler#sample}
		 */
		private Optional<String> sample = Optional.empty();
		/**
		 * See {@link VCFSampler#variantsAmount}
		 */
		private int variantsAmount = 0;
		/**
		 * {@link VCFSampler#afIdentifier}
		 */
		private String afIdentifier = null;
		/**
		 * {@link VCFSampler#acIdentifier}
		 */
		private String acIdentifier = null;
		/**
		 * {@link VCFSampler#anIdentifier}
		 */
		private String anIdentifier = null;
		/**
		 * {@link VCFSampler#filters}
		 */
		private ImmutableSet<IFilter> filters = ImmutableSet.<IFilter> builder().build();
		/**
		 * {@link VCFSampler#intervals}
		 */
		IntervalList intervals = new IntervalList(new SAMFileHeader());
		/**
		 * {@link VCFSampler#deNovoGenerator}
		 */
		DeNovoSampler deNovoSampler;
		/**
		 * The seed to set the {@link Random} number generator for {@link VCFSampler#random}
		 */
		private long seed;
		/**
		 * If <code>true</code> {@link Random} is initialized with the given seed. Otherwise the default {@link Random}
		 * constructor is used.
		 */
		private boolean useSeed = false;

		/**
		 * Builder constructor. Please set at least the {@link #reader} using {@link #vcfReader(VCFFileReader)}.
		 * {@link VCFSampler} can be build running {@link #build()}.
		 */
		public Builder() {
		}

		/**
		 * Create the {@link VCFFileReader} using the path to the file
		 * 
		 * @param path
		 *            The file path to the vcf (bgzip and indexed)
		 * @return The builder with a created {@link #reader}.
		 */
		public Builder file(String path) {
			this.reader = new VCFFileReader(new File(path));
			return this;
		}

		/**
		 * Create the {@link VCFFileReader} using a {@link File}.
		 * 
		 * @param file
		 *            The VCF file (bgzip and indexed)
		 * @return The builder with a created {@link #reader}.
		 */
		public Builder file(File file) {
			this.reader = new VCFFileReader(file);
			return this;
		}

		/**
		 * Set the {@link VCFFileReader}
		 * 
		 * @param reader
		 *            The file of the vcf file
		 * @return The builder with the initialized {link #reader}.
		 */
		public Builder vcfReader(VCFFileReader reader) {
			this.reader = reader;
			return this;
		}

		/**
		 * @param probability
		 *            The probability to select a variant. If set to 1 all variants will be selected.
		 * @return The builder with the initialized {@link #probability}
		 */
		public Builder probability(double probability) {
			this.probability = probability;
			return this;
		}

		/**
		 * @param sample
		 *            If set the genotypes of this samples are selected. if null no sample is selected.
		 * 
		 * @return The builder with the initialized {@link #sample}
		 */
		public Builder sample(String sample) {
			this.sample = Optional.of(sample);
			return this;
		}

		/**
		 * @param afIdentifier
		 *            The identifier of the allele frequency in the Info column.
		 * @return The builder with the initialized {@link #afIdentifier}
		 */
		public Builder afIdentifier(String afIdentifier) {
			this.afIdentifier = afIdentifier;
			return this;
		}

		/**
		 * @param acIdentifier
		 *            The identifier of the alt allele counts in the Info column.
		 * @return The builder with the initialized {@link #acIdentifier}
		 */
		public Builder acIdentifier(String acIdentifier) {
			this.acIdentifier = acIdentifier;
			return this;
		}

		/**
		 * @param anIdentifier
		 *            The identifier of the total number of alleles.
		 * @return The builder with the initialized {@link #anIdentifier}
		 */
		public Builder anIdentifier(String anIdentifier) {
			this.anIdentifier = anIdentifier;
			return this;
		}

		/**
		 * @param variantsAmount
		 *            The maximum number of variants that should be selected by the sampler. If set to smaller than 1
		 *            this option is of and there will be no maximum limit used.
		 * @return The builder with the initialized {@link #variantsAmount}
		 */
		public Builder variantsAmount(int variantsAmount) {
			this.variantsAmount = variantsAmount;
			return this;
		}

		/**
		 * @param intervals
		 *            For each interval the {@link #reader} is queried using the intervals and only such variants are
		 *            parsed
		 * @return The builder with the initialized {@link #intervals}
		 */
		public Builder intervals(IntervalList intervals) {
			this.intervals = intervals.uniqued().sorted();
			return this;
		}

		/**
		 * @param deNovoSampler
		 * @return The builder with the initialized {@link #deNovoSampler}
		 */
		public Builder deNovoGenerator(DeNovoSampler deNovoSampler) {
			this.deNovoSampler = deNovoSampler;
			return this;
		}

		/**
		 * @param filters
		 *            Set of filters that will be used for each {@link VariantContext} read by the
		 *            {@link VCFFileReader#reader}
		 * @return The builder with the initialized {@link #filters}
		 */
		public Builder filters(ImmutableSet<IFilter> filters) {
			this.filters = filters;
			return this;
		}

		/**
		 * @param seed
		 *            The seed to set the {@link Random} number generator for {@link VCFSampler#random}
		 * @return The builder with the initialized {@link #seed}
		 */
		public Builder seed(long seed) {
			this.seed = seed;
			this.useSeed = true;
			return this;
		}

		/**
		 * Build the {@link VCFSampler}. Important: the {@link #reader} have to be initialized using
		 * {@link #file(File)}, {@link #file(String)} or {@link #vcfReader(VCFFileReader)}.
		 * 
		 * @return The created {@link VCFSampler}.
		 */
		public VCFSampler build() {

			// returns null if the reader is not set (not useful)
			if (reader == null)
				throw new RuntimeException("No variants are set for the sampler: The reader cannot be null");

			// do not use it if one of them is not set!
			if (acIdentifier == null || anIdentifier == null) {
				acIdentifier = null;
				anIdentifier = null;
			}

			// check if identifiers ar present

			List<String> ids = Lists.newArrayList(acIdentifier, anIdentifier, afIdentifier);
			ids.removeIf(Objects::isNull);
			for (String id : ids) {
				if (reader.getFileHeader().getInfoHeaderLine(id) == null) {
					throw new InfoIDNotFoundException(id);
				}
			}

			if (variantsAmount < 0)
				variantsAmount = 0;
			else if (variantsAmount > 0) {
				VCFAlternativeAlleleCounter counter = new VCFAlternativeAlleleCounter(reader.iterator(), filters,
						sample);
				setCounts(counter.getCounts());
			}

			Random random;
			if (useSeed)
				random = new Random(seed);
			else
				random = new Random();

			return new VCFSampler(reader, probability, sample, variantsAmount, afIdentifier, acIdentifier, anIdentifier,
					filters, intervals, deNovoSampler, selectAlleles, random);
		}

		private void setCounts(int counts) {
			List<Integer> randomAlleles = new ArrayList<Integer>(counts);
			for (int i = 0; i < counts; i++) {
				randomAlleles.add(i);
			}
			Collections.shuffle(randomAlleles);
			this.selectAlleles = new ArrayList<Integer>(this.variantsAmount);
			for (int i = 0; i < this.variantsAmount; i++) {
				this.selectAlleles.add(randomAlleles.get(i));
			}
			Collections.sort(this.selectAlleles);
		}

	}

	/**
	 * Returns the iterator. If {@link #iterator} is <code>null</code> then the iterator of the {@link #reader} ill be
	 * set. But if intervals are used is generates with {@link #getNextIntervalInterator()} the iterator querying the
	 * next interval.
	 * 
	 * @return The actual VCF file iterator.
	 */
	private CloseableIterator<VariantContext> getIterator() {

		if (this.iterator == null) {
			if (useIntervals())
				this.iterator = getNextIntervalInterator();
			else
				this.iterator = this.reader.iterator();
		}
		return iterator;
	}

	private CloseableIterator<VariantContext> getNextIntervalInterator() {
		Interval interval = nextInterval();
		if (interval != null)
			return this.reader.query(interval.getContig(), interval.getStart(), interval.getEnd());
		else
			return null;

	}

	private Interval nextInterval() {
		Interval output = null;
		if (intervalPosition < getIntervals().size()) {
			output = getIntervals().getIntervals().get(intervalPosition);
			intervalPosition++;
		}
		return output;
	}

	@Override
	public boolean hasNext() {
		return next != null;
	}

	/**
	 * Has next can be true, but next can give back null! But it is important to set the next iterator.
	 * 
	 * @return boolean
	 */
	private boolean checkForNext() {
		// FIXME
		if (useIntervals())
			while (getIterator() != null && !getIterator().hasNext()) {
				// close the actual
				this.iterator.close();
				// get a new one
				this.iterator = getNextIntervalInterator();
			}
		if (getIterator() == null)
			return false;
		return getIterator().hasNext();
	}

	private boolean useIntervals() {
		return !getIntervals().getIntervals().isEmpty();
	}

	@Override
	public VariantContext next() {
		VariantContext actual = next;
		Optional<VariantContext> vc = getNextVariant();
		if (vc.isPresent())
			next = vc.get();
		else
			next = null;
		return actual;
	}

	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}

	private Optional<VariantContext> getNextVariant() {
		Optional<VariantContext> output = Optional.empty();
		while (checkForNext() && !output.isPresent()) {
			// get next line
			VariantContext candidate = getIterator().next();

			// filter
			Optional<VariantContext> optional_candidate = filter(candidate);
			if (!optional_candidate.isPresent())
				continue;
			candidate = optional_candidate.get();

			// get alleles by sampling method
			Map<Integer, Boolean> alleles = useAlleles(candidate);
			if (!alleles.isEmpty()) {
				output = createVariantContextWithGenotype(candidate, alleles);
				if (!output.isPresent())
					continue;
				break;
			}
		}
		return output;
	}

	/**
	 * Applies all filters to the input variant.
	 * 
	 * @param candidate
	 *            Variant to filter
	 * @return The filtered {@link VariantContext}, or <code>null</code>
	 */
	private Optional<VariantContext> filter(VariantContext candidate) {
		Optional<VariantContext> optional_vc = Optional.of(candidate);
		Stream<Optional<VariantContext>> candidateStream = Stream.of(optional_vc);
		for (IFilter iFilter : getFilters()) {
			candidateStream = candidateStream.map(iFilter::filter);
		}
		return candidateStream.collect(MoreCollectors.onlyElement());
	}

	private Optional<VariantContext> createVariantContextWithGenotype(VariantContext candidate,
			Map<Integer, Boolean> alleles) {
		if (sample.isPresent()) {
			Genotype genotype = candidate.getGenotype(sample.get());
			if (genotype.isMixed() || genotype.isNoCall())
				return Optional.empty();
			else if (!genotype.isHomRef())
				return Optional.of(
						new VariantContextBuilder(candidate).genotypes(candidate.getGenotypes(sample.get())).make());
			else
				return Optional.empty();
		} else {
			return Optional.of(new VariantContextBuilder(candidate)
					.genotypes(createGenotype(candidate.getAlleles(), alleles)).make());
		}
	}

	@SuppressWarnings("unused")
	private Collection<Allele> getAlleles(VariantContext candidate, Set<Integer> posOfAllele) {
		Collection<Allele> alleles = new ArrayList<Allele>();
		alleles.add(candidate.getReference());
		for (int i : posOfAllele) {
			alleles.add(candidate.getAlternateAlleles().get(i));
		}
		return alleles;
	}

	private Genotype createGenotype(List<Allele> alleles, Map<Integer, Boolean> use) {
		List<Allele> filteredAlleles = new ArrayList<Allele>();

		int allele = use.keySet().iterator().next() + 1;
		// more the one alternative allele, do not use ref!
		if (use.size() > 1) {
			for (int i : use.keySet()) {
				filteredAlleles.add(alleles.get(i + 1));
			}
		} else if (use.get(allele - 1)) {

			filteredAlleles.add(alleles.get(allele));
			filteredAlleles.add(alleles.get(allele));
		} else {
			filteredAlleles.add(alleles.get(0));
			filteredAlleles.add(alleles.get(allele));

		}

		return GenotypeBuilder.create(DEFAULT_SAMPLE_NAME, filteredAlleles);
	}

	private Map<Integer, Boolean> useAlleles(VariantContext candidate) {
		Map<Integer, Boolean> candidates = new HashMap<Integer, Boolean>();

		if (useAF()) {// AF flag
			Object af = candidate.getCommonInfo().getAttribute(getAFIdentifier());
			if (af instanceof ArrayList<?>) {
				if (((ArrayList<?>) af).get(0) instanceof String) {
					int i = 0;
					for (Object o : (ArrayList<?>) af) {
						addCandidateByHardyWeinberg(candidates, i, Double.parseDouble((String) o));
						i++;
					}
				}
			} else {
				addCandidateByHardyWeinberg(candidates, 0,
						candidate.getCommonInfo().getAttributeAsDouble(getAFIdentifier(), 0.0));
			}
		} else if (useAC()) {
			Object ac = candidate.getCommonInfo().getAttribute(getACIdentifier());
			int an = candidate.getCommonInfo().getAttributeAsInt(getANIdentifier(), 0);
			if (ac instanceof ArrayList<?>) {
				if (((ArrayList<?>) ac).get(0) instanceof String) {
					int i = 0;
					for (Object o : (ArrayList<?>) ac) {
						addCandidateByHardyWeinberg(candidates, i, Double.parseDouble((String) o) / (double) an);
						i++;
					}
				}
			} else {
				addCandidateByHardyWeinberg(candidates, 0,
						(double) candidate.getCommonInfo().getAttributeAsInt(getACIdentifier(), 0) / (double) an);
			}

		} else if (useCounts()) { // variantsAmount > 0
			for (int i = 0; i < candidate.getAlternateAlleles().size(); i++) {
				this.position++;
				if (selectAlleles.contains(position))
					candidates.put(i, nextDouble() <= 0.5);
			}

		} else { // probability
			for (int i = 0; i < candidate.getAlternateAlleles().size(); i++) {
				addCandidateByHardyWeinberg(candidates, 0, getProbability());
			}
		}
		return candidates;

	}

	private boolean useAC() {
		return getACIdentifier() != null && getANIdentifier() != null;
	}

	private void addCandidateByHardyWeinberg(Map<Integer, Boolean> candidates, int i, double af) {
		double[] homHetAF = getHardyWeinbergPrincipleHomhet(af);
		double random = nextDouble();
		if (random <= homHetAF[0])
			candidates.put(i, true);
		else if (random <= homHetAF[1])
			candidates.put(i, false);
	}

	private double[] getHardyWeinbergPrincipleHomhet(double af) {
		return new double[] { Math.pow(1.0 - Math.sqrt(1.0 - af), 2), af };
	}

	private boolean useCounts() {
		return getVariantsAmount() > 0;
	}

	private double nextDouble() {
		return random.nextDouble();
	}

	/**
	 * Getter of the {@link #sample}. Can be null.
	 * 
	 * @return The set sample, or null
	 */
	public Optional<String> getSample() {
		return sample;
	}

	public double getProbability() {
		return probability;
	}

	public VCFHeader getVCFHeader() {
		Set<VCFHeaderLine> set = new LinkedHashSet<VCFHeaderLine>();
		set.addAll(reader.getFileHeader().getMetaDataInInputOrder());
		set.add(new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String, "Genotype"));
		return new VCFHeader(set, getSampleNames());
	}

	/**
	 * Get the sample name stored in {@link #sample} of the new samples.
	 * 
	 * @return The sample name. if no name is set the default name is {@link VCFSampler#DEFAULT_SAMPLE_NAME}.
	 */
	public ImmutableSet<String> getSampleNames() {
		if (!getSample().isPresent())
			return ImmutableSet.<String> builder().add(DEFAULT_SAMPLE_NAME).build();
		return ImmutableSet.<String> builder().add(getSample().get()).build();
	}

	public int getVariantsAmount() {
		return variantsAmount;
	}

	public String getAFIdentifier() {
		return afIdentifier;
	}

	private boolean useAF() {
		return getAFIdentifier() != null;
	}

	public void close() {
		iterator.close();
		reader.close();
	}

	public ImmutableSet<IFilter> getFilters() {
		return filters;
	}

	public String getACIdentifier() {
		return acIdentifier;
	}

	public String getANIdentifier() {
		return anIdentifier;
	}

	public IntervalList getIntervals() {
		return intervals;
	}

}
