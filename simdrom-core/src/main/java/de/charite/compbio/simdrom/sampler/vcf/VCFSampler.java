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
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Lists;
import com.google.common.collect.MoreCollectors;
import com.google.common.math.IntMath;

import de.charite.compbio.simdrom.exception.InfoIDNotFoundException;
import de.charite.compbio.simdrom.filter.IFilter;
import de.charite.compbio.simdrom.filter.RemoveGonosomesFilter;
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
 * The {@link VCFSampler} class uses a VCF file as basis and acts as an iterator
 * over such a file. There are multiple possible options what the
 * {@link VCFSampler} does with the read {@link VariantContext}s.
 * <p>
 * <dl>
 * <dt>Filter</dt>
 * <dd>Variants can be filtered using multiple {@link IFilter}. Completely
 * removed variants where not finally passed into the iterator.</dd>
 * <dt>Sample selection</dt>
 * <dd>If a multi-VCF file is used there is the possibility to select one sample
 * out of it. It can be done randomly by running the VCFRandomSampleSelecter or
 * you can predefine a sample name. The sample name must be present in the VCF
 * header.</dd>
 * <dt>Variant generation using Hardy-Weinberg principle</dt>
 * <dd>If allele frequencies are present in the pares VCF they can be used to
 * generate a completely new artificial {@link Genotype} by using the
 * Hardy-Weinberg principle. Direct allele frequencies can be used of the allele
 * count divided by the number of total alleles.</dd> *
 * <dt>Interval list</dt>
 * <dd>If an {@link IntervalList} is present only variants of that lists will be
 * pares by the {@link VCFFileReader} and these variants are only used by the
 * {@link VCFSampler}.</dd>
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
	private String sampleName;
	/**
	 * The sex of the sampled individual
	 */
	private Sex sex;
	/**
	 * The probability to select a variant. If set to 1 all variants will be
	 * selected.
	 */
	private double probability;
	/**
	 * 
	 */
	private List<Integer> selectAlleles;
	/**
	 * The maximum number of variants that should be selected by the sampler. If
	 * set to smaller than 1 this option is of and there will be no maximum
	 * limit used.
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
	 * The identifier of the het allele counts in the Info column.
	 */
	private String acHetIdentifier;
	/**
	 * The identifier of the hom allele counts in the Info column.
	 */
	private String acHomIdentifier;
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
	 * The next variantContext. If this is null there will be no further variant
	 * and {@link #hasNext()} is false.
	 */
	private VariantContext next;
	/**
	 * Random number generator. Used if probabilities to choose a variant are
	 * smaller than 1.0
	 */
	private Random random;
	/**
	 * De novo generator TODO not used
	 */
	@SuppressWarnings("unused")
	private DeNovoSampler deNovoGenerator;
	/**
	 * Set of filters that will be used for each {@link VariantContext} read by
	 * the {@link #reader}
	 */
	private ImmutableSet<IFilter> filters;
	/**
	 * If list is set (not null) for each interval the {@link #reader} is
	 * queried using the intervals and only such variants are parsed
	 */
	private IntervalList intervals;
	/**
	 * Points to the actual interval in the interval list. If larger than
	 * {@link IntervalList#size()} then variants from all intervals are parsed.
	 */
	private int intervalPosition = 0;

	/**
	 * Default constructor initializing all fields. Can only be used by the
	 * {@link Builder}.
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
	 * @param acHetIdentifier
	 *            {@link #acHetIdentifier}
	 * @param acHomIdentifier
	 *            {@link #acHomIdentifier}
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
	 * @param sampleName
	 *            {@link #sampleName}
	 */
	private VCFSampler(VCFFileReader reader, double probability, Optional<String> sample, int variantsAmount,
			String afIdentifier, String acIdentifier, String acHetIdentifier, String acHomIdentifier,
			String anIdentifier, ImmutableSet<IFilter> filters, IntervalList intervals, DeNovoSampler deNovoSampler,
			List<Integer> selectAlleles, Random random, String sampleName, Sex sex) {
		this.reader = reader;
		this.probability = probability;
		this.sample = sample;
		this.variantsAmount = variantsAmount;
		this.afIdentifier = afIdentifier;
		this.acIdentifier = acIdentifier;
		this.acHetIdentifier = acHetIdentifier;
		this.acHomIdentifier = acHomIdentifier;
		this.anIdentifier = anIdentifier;
		this.filters = filters;
		this.intervals = intervals;
		this.selectAlleles = selectAlleles;
		this.deNovoGenerator = deNovoSampler;
		this.random = random;
		this.sampleName = sampleName;
		this.sex = sex;

		next();
	}

	/**
	 * 
	 * Builder class of the {@link VCFSampler}. Checks if everything is set
	 * correctly before building. The reader have to be set by using
	 * {@link #vcfReader(VCFFileReader)}.
	 * 
	 * @author <a href="mailto:max.schubach@charite.de">Max Schubach</a>
	 *
	 */
	public static final class Builder {

		/**
		 * See {@link VCFSampler#sex}
		 */
		private Optional<Sex> sex = Optional.empty();

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
		 * {@link VCFSampler#acHetIdentifier}
		 */
		private String acHetIdentifier = null;
		/**
		 * {@link VCFSampler#acHomIdentifier}
		 */
		private String acHomIdentifier = null;
		/**
		 * {@link VCFSampler#anIdentifier}
		 */
		private String anIdentifier = null;
		/**
		 * {@link VCFSampler#filters}
		 */
		private ImmutableSet<IFilter> filters = ImmutableSet.<IFilter>builder().build();
		/**
		 * {@link VCFSampler#intervals}
		 */
		IntervalList intervals = new IntervalList(new SAMFileHeader());
		/**
		 * {@link VCFSampler#deNovoGenerator}
		 */
		DeNovoSampler deNovoSampler;
		/**
		 * The seed to set the {@link Random} number generator for
		 * {@link VCFSampler#random}
		 */
		private long seed;
		/**
		 * If <code>true</code> {@link Random} is initialized with the given
		 * seed. Otherwise the default {@link Random} constructor is used.
		 */
		private boolean useSeed = false;

		private String sampleName = "Sampled";

		/**
		 * Builder constructor. Please set at least the {@link #reader} using
		 * {@link #vcfReader(VCFFileReader)}. {@link VCFSampler} can be build
		 * running {@link #build()}.
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
		 *            The probability to select a variant. If set to 1 all
		 *            variants will be selected.
		 * @return The builder with the initialized {@link #probability}
		 */
		public Builder probability(double probability) {
			this.probability = probability;
			return this;
		}

		/**
		 * @param sample
		 *            If set the genotypes of this samples are selected. if null
		 *            no sample is selected.
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
		 *            The identifier of the alt allele counts in the Info
		 *            column.
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
		 * @param acHetIdentifier
		 *            The identifier of the heterozygouse number of alleles.
		 * @return The builder with the initialized {@link #acHetIdentifier}
		 */
		public Builder acHetIdentifier(String acHetIdentifier) {
			this.acHetIdentifier = acHetIdentifier;
			return this;
		}

		/**
		 * @param acHomIdentifier
		 *            The identifier of the homozygouse number of alleles.
		 * @return The builder with the initialized {@link #acHomIdentifier}
		 */
		public Builder acHomIdentifier(String acHomIdentifier) {
			this.acHomIdentifier = acHomIdentifier;
			return this;
		}

		/**
		 * @param variantsAmount
		 *            The maximum number of variants that should be selected by
		 *            the sampler. If set to smaller than 1 this option is of
		 *            and there will be no maximum limit used.
		 * @return The builder with the initialized {@link #variantsAmount}
		 */
		public Builder variantsAmount(int variantsAmount) {
			this.variantsAmount = variantsAmount;
			return this;
		}

		/**
		 * @param sex
		 *            The sex of the individual that should be sampled. Only
		 *            important if Allele frequencies or allele counts are used.
		 *            For male X-hom alt and Y are sampled. For females only X.
		 *            If set to {@links Sex#UNKNOWN} no gonosomes will be
		 *            sampled!
		 * @return The builder with the initialized {@link #sex}
		 */
		public Builder sex(Sex sex) {
			this.sex = Optional.of(sex);
			return this;
		}

		/**
		 * @param intervals
		 *            For each interval the {@link #reader} is queried using the
		 *            intervals and only such variants are parsed
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
		 *            Set of filters that will be used for each
		 *            {@link VariantContext} read by the
		 *            {@link VCFFileReader#reader}
		 * @return The builder with the initialized {@link #filters}
		 */
		public Builder filters(ImmutableSet<IFilter> filters) {
			this.filters = filters;
			return this;
		}

		/**
		 * @param seed
		 *            The seed to set the {@link Random} number generator for
		 *            {@link VCFSampler#random}
		 * @return The builder with the initialized {@link #seed}
		 */
		public Builder seed(long seed) {
			this.seed = seed;
			this.useSeed = true;
			return this;
		}

		/**
		 * Build the {@link VCFSampler}. Important: the {@link #reader} have to
		 * be initialized using {@link #file(File)}, {@link #file(String)} or
		 * {@link #vcfReader(VCFFileReader)}.
		 * 
		 * @return The created {@link VCFSampler}.
		 */
		public VCFSampler build() {

			// returns null if the reader is not set (not useful)
			if (reader == null)
				throw new RuntimeException("No variants are set for the sampler: The reader cannot be null");

			// do not use it if one of them is not set!
			if (anIdentifier != null) {
				if (acIdentifier != null) {
					acHetIdentifier = null;
					acHomIdentifier = null;
				} else if (acHetIdentifier != null && acHomIdentifier != null)
					acIdentifier = null;
				else {
					acIdentifier = null;
					anIdentifier = null;
					acHetIdentifier = null;
					acHomIdentifier = null;
				}
			} else {
				acIdentifier = null;
				anIdentifier = null;
				acHetIdentifier = null;
				acHomIdentifier = null;
			}

			// check if identifiers are present
			List<String> ids = Lists.newArrayList(acIdentifier, anIdentifier, afIdentifier, acHetIdentifier,
					acHomIdentifier);
			ids.removeIf(Objects::isNull);
			for (String id : ids) {
				if (reader.getFileHeader().getInfoHeaderLine(id) == null) {
					throw new InfoIDNotFoundException(id);
				}
			}

			Random random;
			if (useSeed)
				random = new Random(seed);
			else
				random = new Random();

			// Random sex if not set
			if (!sex.isPresent())
				sex = Optional.of(random.nextBoolean() ? Sex.MALE : Sex.FEMALE);

			// set sex to unknown if no sample is selected and remove then all
			// gonosomes
			if (sample.isPresent()) {
				sex = Optional.of(Sex.UNKNOWN);
			}

			// remove gonosomes if sex is none
			if (sex.get() == Sex.NONE)
				filters = ImmutableSet.<IFilter>builder().add(new RemoveGonosomesFilter()).addAll(filters).build();

			if (variantsAmount < 0)
				variantsAmount = 0;
			else if (variantsAmount > 0) {
				VCFAlternativeAlleleCounter counter = new VCFAlternativeAlleleCounter(reader.iterator(), filters,
						sample);
				setCounts(counter.getCounts(), random);
			}

			return new VCFSampler(reader, probability, sample, variantsAmount, afIdentifier, acIdentifier,
					acHetIdentifier, acHomIdentifier, anIdentifier, filters, intervals, deNovoSampler, selectAlleles,
					random, sampleName, sex.get());
		}

		private void setCounts(int counts, Random random) {
			List<Integer> randomAlleles = new ArrayList<Integer>(counts);
			for (int i = 0; i < counts; i++) {
				randomAlleles.add(i);
			}
			Collections.shuffle(randomAlleles, random);
			this.selectAlleles = new ArrayList<Integer>(this.variantsAmount);
			for (int i = 0; i < this.variantsAmount; i++) {
				this.selectAlleles.add(randomAlleles.get(i));
			}
			Collections.sort(this.selectAlleles);
		}

		public void sampleName(String sampleName) {
			this.sampleName = sampleName;

		}

	}

	/**
	 * Returns the iterator. If {@link #iterator} is <code>null</code> then the
	 * iterator of the {@link #reader} will be set. But if intervals are used is
	 * generates with {@link #getNextIntervalInterator()} the iterator querying
	 * the next interval.
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
	 * Has next can be true, but next can give back null! But it is important to
	 * set the next iterator.
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
			Optional<int[]> alleles = useAlleles(candidate);
			if (alleles.isPresent()) {
				output = createVariantContextWithGenotype(candidate, alleles.get());
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

	private Optional<VariantContext> createVariantContextWithGenotype(VariantContext candidate, int[] alleles) {
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
			boolean onlyRef = true;
			for (int i : alleles) {
				if (i != 0) {
					onlyRef = false;
					break;
				}
			}
			if (onlyRef)
				return  Optional.empty();
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

	private Genotype createGenotype(List<Allele> alleles, int[] use) {
		List<Allele> filteredAlleles = new ArrayList<Allele>();
		for (int allele : use) {
			filteredAlleles.add(alleles.get(allele));
		}

		return GenotypeBuilder.create(sampleName, filteredAlleles);
	}

	private Optional<int[]> useAlleles(VariantContext candidate) {
		Optional<int[]> candidates = Optional.empty();
		double[] afs;
		if (useAF()) {// AF flag
			Object af = candidate.getCommonInfo().getAttribute(getAFIdentifier());
			if (af instanceof ArrayList<?>) {
				ArrayList<Double> afAltList = new ArrayList<>();
				if (((ArrayList<?>) af).get(0) instanceof String) {
					for (Object o : (ArrayList<?>) af) {
						afAltList.add(Double.parseDouble((String) o));
					}
				}
				afs = createAFsFromAltAFs(afAltList);
			} else {
				double altAF = candidate.getCommonInfo().getAttributeAsDouble(getAFIdentifier(), 0.0);
				afs = new double[] { 1.0 - altAF, altAF };
			}
			candidates = getCandidateByHardyWeinberg(candidate.getContig(), afs);
		} else if (useAC()) {
			Object ac = candidate.getCommonInfo().getAttribute(getACIdentifier());
			int an = candidate.getCommonInfo().getAttributeAsInt(getANIdentifier(), 0);
			if (ac instanceof ArrayList<?>) {
				ArrayList<Double> afAltList = new ArrayList<>();
				if (((ArrayList<?>) ac).get(0) instanceof String) {
					for (Object o : (ArrayList<?>) ac) {
						afAltList.add(Double.parseDouble((String) o) / (double) an);
					}
				}
				afs = createAFsFromAltAFs(afAltList);
			} else {
				double altAF = (double) candidate.getCommonInfo().getAttributeAsInt(getACIdentifier(), 0) / (double) an;
				afs = new double[] { 1.0 - altAF, altAF };
			}
		} else if (useACHetHom()) {
			Object acHet = candidate.getCommonInfo().getAttribute(getACHetIdentifier());
			Object acHom = candidate.getCommonInfo().getAttribute(getACHomIdentifier());
			int an = candidate.getCommonInfo().getAttributeAsInt(getANIdentifier(), 0);
			double anIndividuals = an / 2.0;
			if (acHet instanceof ArrayList<?> && acHom instanceof ArrayList<?>) {
				ArrayList<Double> afAltHomList = new ArrayList<>();
				ArrayList<Double> afAltHetList = new ArrayList<>();
				if (((ArrayList<?>) acHet).get(0) instanceof String
						&& ((ArrayList<?>) acHom).get(0) instanceof String) {
					for (Object oHom : (ArrayList<?>) acHom) {
						afAltHomList.add(Double.parseDouble((String) oHom) / anIndividuals);
					}
					for (Object oHet : (ArrayList<?>) acHet) {
						afAltHetList.add(Double.parseDouble((String) oHet) / anIndividuals);
					}
				}
				double sum = 0;
				double[] afsHom = new double[afAltHomList.size() + 1];
				double[] afsHet = new double[IntMath.binomial(afAltHomList.size(), 2)];
				int pos = 0;
				for (int i = 0; i < candidate.getAlleles().size(); i++) {
					for (int j = i; i < candidate.getAlleles().size(); j++) {
						if (i == j) {
							sum += afAltHomList.get(i);
							afsHom[i + 1] = afAltHomList.get(i);
						} else {
							sum += afAltHetList.get(pos);
							afsHet[pos] = afAltHetList.get(pos);
							pos++;
						}
					}
				}
				afsHom[0] = Math.max(0.0, 1.0 - sum);
				candidates = getCandidate(new double[][] { afsHom, afsHet });
			} else {
				double altHomAF = (double) candidate.getCommonInfo().getAttributeAsInt(getACHomIdentifier(), 0)
						/ anIndividuals;
				double altHetAF = (double) candidate.getCommonInfo().getAttributeAsInt(getACHetIdentifier(), 0)
						/ anIndividuals;
				candidates = getCandidate(new double[][] { new double[] { 1.0 - altHomAF - altHetAF, altHomAF },
						new double[] { altHetAF } });
			}
		} else if (useCounts()) { // variantsAmount > 0
			for (int i = 0; i < candidate.getAlternateAlleles().size(); i++) {
				this.position++;
				if (selectAlleles.contains(position)) {
					if (nextDouble() <= 0.25)
						candidates = Optional.of(new int[] { i + 1, i + 1 });
					else
						candidates = Optional.of(new int[] { 0, i + 1 });
					break;
				}
			}

		} else { // probability
			afs = new double[candidate.getAlleles().size()];
			afs[0] = 1.0 - getProbability();
			int altSize = candidate.getAlternateAlleles().size();
			for (int i = 1; i < afs.length; i++) {
				afs[i] = getProbability() / (double) altSize;
			}
			candidates = getCandidateByHardyWeinberg(candidate.getContig(), afs);
		}
		return candidates;

	}

	private double[] createAFsFromAltAFs(ArrayList<Double> afAltList) {
		double[] afs;
		double afRef = 1.0;
		for (Double afAlt : afAltList) {
			afRef -= afAlt;
		}
		afs = new double[afAltList.size() + 1];
		afs[0] = Math.max(0, afRef);
		for (int i = 0; i < afAltList.size(); i++) {
			afs[i + 1] = afAltList.get(i);
		}
		return afs;
	}

	/**
	 * @param contig
	 * @return
	 */
	private boolean onX(String contig) {
		return contig.matches("chr23|chrX|Chr23|ChrX|X|x|23|chrx|Chrx");
	}

	/**
	 * @param contig
	 * @return
	 */
	private boolean onY(String contig) {
		return contig.matches("chr24|chrY|Chr24|ChrY|Y|y|24|chry|Chry");
	}

	/**
	 * @param candidates
	 * @param i
	 * @param candidate
	 */
	private void addPossibleCandidateAllele(Map<Integer, Boolean> candidates, int i,
			Optional<Boolean> candidateAllele) {
		if (candidateAllele.isPresent())
			candidates.put(i, candidateAllele.get());
	}

	private Map<Integer, Boolean> selectedMaxTwoAlleles(Map<Integer, Boolean> candidates) {
		Map<Integer, Boolean> output = new HashMap<>();
		int count = 0;
		List<Integer> integers = Lists.newArrayList(candidates.keySet());
		Collections.shuffle(integers, random);

		for (Integer genotype : integers) {
			if (count == 2)
				return output;
			int alleles = (candidates.get(genotype) ? 2 : 1);
			if (count + alleles > 2) {
				count += 1;
				output.put(genotype, false);
			} else {
				count += alleles;
				output.put(genotype, candidates.get(genotype));
			}

		}
		return output;
	}

	private int countAlleles(Map<Integer, Boolean> candidates) {
		int i = 0;
		for (Boolean hom : candidates.values()) {
			if (hom)
				i += 2;
			else
				i += 1;
		}
		return i;
	}

	private boolean useAC() {
		return getACIdentifier() != null && getANIdentifier() != null;
	}

	private boolean useACHetHom() {
		return getACHetIdentifier() != null && getACHomIdentifier() != null && getANIdentifier() != null;
	}

	private Optional<int[]> getCandidateByHardyWeinberg(String contig, double[] af) {
		boolean onX = onX(contig);
		boolean onY = onY(contig);
		double[][] homHetProbs;
		if (onX | onY) {

			switch (sex) {
			case MALE:
				if (onX)
					homHetProbs = getHardyWeinbergPrincipleMaleXHom(af);
				else
					homHetProbs = getHardyWeinbergPrincipleMaleYHom(af);
				break;

			case FEMALE:
				if (onX)
					homHetProbs = getHardyWeinbergPrincipleFemaleXHomhet(af);
				else
					homHetProbs = new double[][] { new double[0], new double[0] };
				break;

			case UNKNOWN:
				if (onX)
					homHetProbs = (random.nextBoolean() ? getHardyWeinbergPrincipleMaleXHom(af)
							: getHardyWeinbergPrincipleFemaleXHomhet(af));
				else
					homHetProbs = (random.nextBoolean() ? getHardyWeinbergPrincipleMaleYHom(af)
							: new double[][] { new double[0], new double[0] });
				break;
			case NONE:

				return Optional.empty();
			}
		}

		homHetProbs = getHardyWeinbergPrincipleHomhet(af);

		return getCandidate(homHetProbs);
	}

	/**
	 * @param af
	 * @return
	 */
	private double[][] getHardyWeinbergPrincipleFemaleXHomhet(double[] af) {
		double[] homs = new double[af.length];
		double[] hets = new double[IntMath.binomial(af.length, 2)];
		for (int i = 0; i < af.length; i++) {
			homs[i] = Math.pow(af[i], 2);
		}
		for (int i = 0; i < af.length; i++) {
			for (int j = i + 1; j < af.length; j++) {
				hets[i] = 2.0 * af[i] * af[j];
			}
		}
		return new double[][] { homs, hets };
	}

	/**
	 * @param af
	 * @return
	 */
	private double[][] getHardyWeinbergPrincipleMaleYHom(double[] af) {
		double[] homs = new double[af.length];
		double[] hets = new double[IntMath.binomial(af.length, 2)];
		for (int i = 0; i < af.length; i++) {
			homs[i] = Math.pow(af[i], 2);
		}
		for (int i = 0; i < af.length; i++) {
			for (int j = i + 1; j < af.length; j++) {
				hets[i] = 2.0 * af[i] * af[j];
			}
		}
		return new double[][] { homs, hets };
	}

	/**
	 * @param af
	 * @return
	 */
	private double[][] getHardyWeinbergPrincipleMaleXHom(double[] af) {
		double[] homs = new double[af.length];
		double[] hets = new double[IntMath.binomial(af.length, 2)];
		for (int i = 0; i < af.length; i++) {
			homs[i] = Math.pow(af[i], 2);
		}
		for (int i = 0; i < af.length; i++) {
			for (int j = i + 1; j < af.length; j++) {
				hets[i] = 2.0 * af[i] * af[j];
			}
		}
		return new double[][] { homs, hets };
	}

	private Optional<int[]> getCandidate(double[][] homHetProbs) {
		double random = nextDouble();
		double sum = 0.0;
		for (int i = 0; i < homHetProbs[0].length; i++) {
			sum += homHetProbs[0][i];
			if (random <= sum)
				return Optional.of(new int[] { i, i });
		}
		int pos = 0;
		for (int i = 0; i < homHetProbs[0].length; i++) {
			for (int j = i + 1; j < homHetProbs[0].length; j++) {
				sum += homHetProbs[0][pos];
				if (random <= sum)
					return Optional.of(new int[] { i, j });
				pos++;

			}
		}
		return Optional.empty();
	}

	private double[][] getHardyWeinbergPrincipleHomhet(double[] af) {
		double[] homs = new double[af.length];
		double[] hets = new double[IntMath.binomial(af.length, 2)];
		for (int i = 0; i < af.length; i++) {
			homs[i] = Math.pow(af[i], 2);
		}
		for (int i = 0; i < af.length; i++) {
			for (int j = i + 1; j < af.length; j++) {
				hets[i] = 2.0 * af[i] * af[j];
			}
		}
		return new double[][] { homs, hets };
	}

	private boolean useCounts() {
		return getVariantsAmount() > 0;
	}

	private double nextDouble() {
		return random.nextDouble();
	}

	/**
	 * Getter of the {@link #sample}.
	 * 
	 * @return The set sample
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
	 * @return The sample name. if no name is set the default name is
	 *         {@link VCFSampler#sampleName}.
	 */
	public ImmutableSet<String> getSampleNames() {
		if (!getSample().isPresent())
			return ImmutableSet.<String>builder().add(sampleName).build();
		return ImmutableSet.<String>builder().add(getSample().get()).build();
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

	public String getACHetIdentifier() {
		return acHetIdentifier;
	}

	public String getACHomIdentifier() {
		return acHomIdentifier;
	}

	public String getANIdentifier() {
		return anIdentifier;
	}

	public IntervalList getIntervals() {
		return intervals;
	}

	public String getSampleName() {
		return sampleName;
	}

}
